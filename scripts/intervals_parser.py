from datetime import date
import gcsfs
import pandas as pd

import hydra

from pyliftover import LiftOver

from pyspark.sql import dataframe
import pyspark.sql
import pyspark.sql.types as t
import pyspark.sql.functions as f
from pyspark.sql import SparkSession


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
        self.liftover_udf = f.udf(lambda chrom, pos: self.lo.convert_coordinate(chrom, pos), t.ArrayType(t.ArrayType(t.StringType())))

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
                f.col('mapped_start').isNotNull() & f.col('mapped_end').isNotNull()

                # Drop rows where the start is larger than the end:
                & (f.col('mapped_' + end_col) >= f.col('mapped_' + start_col))

                # Drop rows where the difference of the length of the regions are larger than the threshold:
                & (
                    f.abs(
                        (f.col('mapped_' + end_col) - f.col('mapped_' + start_col)) -
                        (f.col('mapped_' + end_col) - f.col('mapped_' + start_col))
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
            .withColumn('mapped', self.liftover_udf(f.col(chrom_name), f.col(pos_name)))
            .filter((f.col('mapped').isNotNull()) & (f.size(f.col('mapped')) == 1))

            # Extracting mapped corrdinates:
            .withColumn('mapped_' + chrom_name, f.col('mapped')[0][0])
            .withColumn('mapped_' + pos_name, f.col('mapped')[0][1])

            # Drop rows that mapped to the other chromosomes:
            .filter(f.col('mapped_' + chrom_name) == f.concat(f.lit('chr'),f.col(chrom_name)))

            # Dropping unused coluns:
            .drop('mapped', 'mapped_' + chrom_name)
            .persist()
        )

        return mapped

def parse_anderson(cfg, gene_index, lift, proximity_limit):
    """
    Parse the anderson file and return a dataframe with the intervals.

    :param anderson_file: Path to the anderson file (.bed).
    :return: Spark Dataframe with the intervals.
    """

    # Read the anderson file:
    anderson_df = (
        SparkSession.getActiveSession().createDataFrame(pd.read_csv(cfg.data_file, sep='\t', header=0, low_memory=False, skiprows=1))
        # Parsing score column and casting as float:
        .withColumn('score', f.col('score').cast('float') / f.lit(1000))

        # Parsing the 'name' column:
        .withColumn('parsedName', f.split(f.col('name'), ';'))
        .withColumn('gene_symbol', f.col('parsedName')[2])
        .withColumn('location', f.col('parsedName')[0])
        .withColumn('chrom', f.regexp_replace(f.split(f.col('location'), ':|-')[0], 'chr', ''))
        .withColumn('start', f.split(f.col('location'), ':|-')[1].cast(t.IntegerType()))
        .withColumn('end', f.split(f.col('location'), ':|-')[2].cast(t.IntegerType()))

        # Select relevant columns:
        .select('chrom', 'start', 'end', 'gene_symbol', 'score')

        # Drop rows with non-canonical chromosomes:
        .filter(f.col('chrom').isin([str(x) for x in range(1, 23)] + ['X', 'Y', 'MT']))

        # For each region/gene, keep only one row with the highest score:
        .groupBy('chrom', 'start', 'end', 'gene_symbol')
        .agg(f.max('score').alias('score'))

        .orderBy('chrom', 'start')
        .persist()
    )

    # Selecting relevant columns from gene set:
    genes = gene_index.withColumnRenamed('gene_name', 'gene_symbol').select('gene_symbol', 'chr', 'gene_id', 'TSS')

    return (
        # Lift over the intervals:
        lift.convert_intervals(anderson_df, 'chrom', 'start', 'end')
        .drop('start', 'end')
        .withColumnRenamed('mapped_start', 'start')
        .withColumnRenamed('mapped_end', 'end')
        .distinct()

        # Joining with the gene index (unfortunately we are losing a bunch of genes here due to old symols):
        .join(genes, on='gene_symbol', how='left')
        .filter(
            # Drop rows where the gene is not on the same chromosome
            (f.col('chrom') == f.regexp_replace(f.col('chr'), 'chr', ''))
            # Drop rows where the TSS is far from the start of the region
            & (f.abs((f.col('start') - f.col('end') / 2) - f.col('TSS')) <= proximity_limit)
        )

        # Adding constant values:
        .withColumn('dataset_name', f.lit(cfg.dataset_name))
        .withColumn('data_type', f.lit(cfg.data_type))
        .withColumn('experiment_type', f.lit(cfg.experiment_type))
        .withColumn('pmid', f.lit(cfg.pubmed_id))
        .withColumn('bio_feature', f.lit(cfg.bio_feature))

        # Select relevant columns:
        .select(
            'chrom', 'start', 'end', 'gene_id', 'score', 'dataset_name', 'data_type', 'experiment_type', 'pmid', 'bio_feature'
        )
        .persist()
    )

@hydra.main(config_path='../configs', config_name='config')
def main(cfg):

    spark = (
        pyspark.sql.SparkSession
        .builder
        .master("local[*]")
        .getOrCreate()
    )

    chain_file = cfg.intervals.liftover_chain_file
    max_difference = cfg.intervals.max_lenght_difference
    proximity_limit = 2 * float(cfg.intervals.proximity_limit)

    # Open and process gene file:
    gene_index = spark.read.parquet(cfg.intervals.gene_index).persist()

    # Initialize liftover object:
    lift = LiftOverSpark(chain_file, max_difference)

    # Parsing the anderson file:
    anderson_df = parse_anderson(cfg.intervals.anderson, gene_index, lift, proximity_limit)

    # Further parsers will come here...

    # Saving data:
    version = date.today().strftime("%y%m%d")
    anderson_df.write.parquet(cfg.intervals.output + f'/interval_{version}')


if __name__ == '__main__':

    # Calling parsers:
    main()
