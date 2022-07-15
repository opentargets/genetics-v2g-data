import gcsfs

from pyliftover import LiftOver

import pyspark.sql.types as T
import pyspark.sql.functions as F
from pyspark.sql import dataframe

class LiftOverSpark:
    """
    LiftOver class for mapping genomic coordinates to an other genome build.

    The input is a Spark DataFrame with a chromosome and position column. This classs can
    also map regions, if a start and end positions are provided.

    **Logic**:

    - The mapping is dropped if the mapped chromosome is not on the same as the source.
    - The mapping is dropped if the mapping is ambiguous (more than one mapping is available).
    - If regions are provided, the mapping is dropped if the new region is reversed (mapped_start > mapped_end).
    - If regions are provided, the mapping is dropped if the difference of the lenght of the mapped region and original is larger than a threshold.
    - When lifting over intervals, only unique coordinates are lifted, they joined back to the original dataframe.
    """
    def __init__(self, chain_file: str, max_difference: int = None) -> None:
        """

        :param chain_file: Path to the chain file. Local or google bucket. Chainfile is not gzipped!
        :param max_difference: Maximum difference between the length of the mapped region and the original region.
        """

        self.chain_file = chain_file

        # Initializing liftover object by opening the chain file:
        if chain_file.startswith('gs://'):
            with gcsfs.GCSFileSystem().open(chain_file) as chain_file_object:
                self.lo = LiftOver(chain_file_object)
        else:
            self.lo = LiftOver(chain_file)

        self.max_difference = max_difference

        # UDF to do map genomic coordinates to liftover coordinates:
        self.liftover_udf = F.udf(lambda chrom, pos: self.lo.convert_coordinate(chrom, pos), T.ArrayType(T.ArrayType(T.StringType())))

    def convert_intervals(self, df: dataframe, chrom_col: str, start_col: str, end_col: str) -> dataframe:
        """
        Convert genomic intervals to liftover coordinates

        :param df: spark Dataframe with chromosome, start and end columns.
        :param chrom_col: Name of the chromosome column.
        :param start_col: Name of the start column.
        :param end_col: Name of the end column.

        :return: filtered Spark Dataframe with the mapped start and end coordinates.
        """

        # Lift over start coordinates:
        start_df = df.select(chrom_col, start_col).distinct()
        start_df = self.convert_coordinates(start_df, chrom_col, start_col).withColumnRenamed('mapped_pos', 'mapped_' + start_col)

        # Lift over end coordinates:
        end_df = df.select(chrom_col, end_col).distinct()
        end_df = self.convert_coordinates(end_df, chrom_col, end_col).withColumnRenamed('mapped_pos', 'mapped_' + end_col)

        # Join dataframe with mappings:
        return (
            df
            .join(start_df, on=[chrom_col, start_col], how='left')
            .join(end_df, on=[chrom_col, end_col], how='left')

            # Select only rows where the start is smaller than the end:
            .filter(
                # Drop rows with no mappings:
                F.col('mapped_' + start_col).isNotNull() & F.col('mapped_' + end_col).isNotNull()

                # Drop rows where the start is larger than the end:
                & (f.col('mapped_' + end_col) >= F.col('mapped_' + start_col))

                # Drop rows where the difference of the length of the regions are larger than the threshold:
                & (
                    F.abs(
                        (F.col('mapped_' + end_col) - F.col('mapped_' + start_col)) -
                        (F.col('mapped_' + end_col) - F.col('mapped_' + start_col))
                    ) <= self.max_difference
                )
            )
            .persist()
        )

    def convert_coordinates(self, df: dataframe, chrom_name: str, pos_name: str) -> list:
        """
        Converts genomic coordinates to coordinates on an other build

        :param df: Spark Dataframe with chromosome and position columns.
        :param chrom_name: Name of the chromosome column.
        :param pos_name: Name of the position column.

        :return: Spark Dataframe with the mapped position column.
        """
        mapped =  (
            df
            .withColumn('mapped', self.liftover_udf(F.col(chrom_name), F.col(pos_name)))
            .filter((F.col('mapped').isNotNull()) & (F.size(F.col('mapped')) == 1))

            # Extracting mapped corrdinates:
            .withColumn('mapped_' + chrom_name, F.col('mapped')[0][0])
            .withColumn('mapped_' + pos_name, F.col('mapped')[0][1])

            # Drop rows that mapped to the other chromosomes:
            .filter(F.col('mapped_' + chrom_name) == F.concat(F.lit('chr'), F.col(chrom_name)))

            # Dropping unused columns:
            .drop('mapped', 'mapped_' + chrom_name)
            .persist()
        )

        return mapped