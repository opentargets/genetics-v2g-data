import argparse
from functools import reduce
import logging
import os

import gcsfs
import pandas as pd

import pyspark
import pyspark.sql.types as T
import pyspark.sql.functions as F
from pyspark.sql import dataframe, SparkSession

from modules.Liftover import LiftOverSpark


class parser_javierre:
    """
    Parser Javierre 2016 dataset
    """

    # Constants:
    DATASET_NAME = 'javierre2016'
    DATA_TYPE = 'interval'
    EXPERIMENT_TYPE = 'pchic'
    PMID = '27863249'
    TWOSIDED_THRESHOLD = 2.45e6

    def __init__(self,
                 javierre_parquet: str,
                 gene_index: dataframe,
                 lift: LiftOverSpark) -> None:

        # Read gene index:
        genes = self.prepare_gene_index(gene_index)

        # Read Javierre data:
        javierre_raw = (
            SparkSession.getActiveSession().read.parquet(javierre_parquet)

            # SPlitting name column into chromosome, start, end, and score:
            .withColumn('name_split', F.split(F.col('name'), r':|-|,'))
            .withColumn('name_chr', F.regexp_replace(F.col('name_split')[0], 'chr', '').cast(T.StringType()))
            .withColumn('name_start', F.col('name_split')[1].cast(T.IntegerType()))
            .withColumn('name_end', F.col('name_split')[2].cast(T.IntegerType()))
            .withColumn('name_score', F.col('name_split')[3].cast(T.FloatType()))

            # Cleaning up chromosome:
            .withColumn('chrom', F.regexp_replace(F.col('chrom'), 'chr', '').cast(T.StringType()))
            .drop('name_split', 'name', 'annotation')

            # Keep canonical chromosomes and consistent chromosomes:
            .filter(
                (F.col('chrom') == F.col('name_chr')) & # 6_473_184 -> 6_428_824
                F.col('name_chr').isin([f'{x}' for x in range(1, 23)] + ['X', 'Y', 'MT'])
            )
        )

        # Lifting over intervals:
        javierre_remapped = (
            javierre_raw
            # Mapping interval 1:
            .transform(lambda df: lift.convert_intervals(df, 'chrom', 'start', 'end', filter=False))
            .drop('start', 'end')
            .withColumnRenamed('mapped_chrom', 'chrom')
            .withColumnRenamed('mapped_start', 'start')
            .withColumnRenamed('mapped_end', 'end')

            # Mapping interval 2:
            .transform(lambda df: lift.convert_intervals(df, 'name_chr', 'name_start', 'name_end', filter=False))
            .drop('name_start', 'name_end')
            .withColumnRenamed('mapped_name_chr', 'name_chr')
            .withColumnRenamed('mapped_name_start', 'name_start')
            .withColumnRenamed('mapped_name_end', 'name_end')

            .persist()
        )

        # Mapping intervals to genes:
        self.javierre_intervals = (
            javierre_remapped
            .join(genes, on=[
                (genes.gene_chr == javierre_remapped.chrom) &
                ((genes.gene_start <= javierre_remapped.end) & (genes.gene_start >= javierre_remapped.start)) |
                ((genes.gene_end <= javierre_remapped.end) & (genes.gene_end >= javierre_remapped.start))
            ])
            .filter(
                # Drop rows where the TSS is far from the start of the region
                F.abs((F.col('start') + F.col('end')) / 2 - F.col('TSS')) <= self.TWOSIDED_THRESHOLD
            )

            # For each gene, keep only the highest scoring interval:
            .groupBy('name_chr', 'name_start', 'name_end', 'gene_id', 'bio_feature')
            .agg(F.max(F.col('name_score')).alias('score'))

            # Create the output:
            .select(
                F.col('name_chr').alias('chr'),
                F.col('name_start').alias('start'),
                F.col('name_end').alias('end'),
                F.col('score'),
                F.col('gene_id').alias('gene_id'),
                F.col('bio_feature').alias('cell_type'),
                F.lit(None).alias('bio_feature'),
                F.lit(self.DATASET_NAME).alias('dataset'),
                F.lit(self.DATA_TYPE).alias('data_type'),
                F.lit(self.EXPERIMENT_TYPE).alias('experiment_type'),
                F.lit(self.PMID).alias('pmid')
            )
            .persist()
        )

    def get_intervals_data(self) -> dataframe:
        return self.javierre_intervals

    def qc_intervals(self) -> None:
        """
        Perform QC on the anderson intervals.
        """

        # Get numbers:
        logging.info(f'Size of Javierre data: {self.javierre_intervals.count()}')
        logging.info(f'Number of unique intervals: {self.javierre_intervals.select("start", "end").distinct().count()}')
        logging.info(f'Number genes in the Javierre dataset: {self.javierre_intervals.select("gene_id").distinct().count()}')

    def save_parquet(self, output_file: str) -> None:
        self.javierre_intervals.write.mode('overwrite').parquet(output_file)

    @staticmethod
    def prepare_gene_index(gene_index: dataframe) -> dataframe:
        '''Pre-processind the gene dataset
        - selecting and renaming relevant columns
        - remove 'chr' from chromosome column

        :param gene_index: Path to the gene parquet file
        :return: Spark Dataframe
        '''
        # Reading gene annotations:
        return (
            gene_index
            .select(
                F.regexp_replace(F.col('chr'), 'chr', '').alias('chr'),
                F.col('start').alias('gene_start'),
                F.col('end').alias('gene_end'),
                'gene_id', 'TSS'
            )
            .withColumn('chr', F.regexp_replace(F.col('chr'), 'chr', ''))
            .withColumnRenamed('chr', 'gene_chr')
            .withColumnRenamed('start', 'gene_start')
            .withColumnRenamed('end', 'gene_end')
            .persist()
        )


def main(javierre_data_file: str, gene_index_file: str, chain_file: str, output_file: str) -> None:

    spark = (
        pyspark.sql.SparkSession
        .builder
        .master("local[*]")
        .getOrCreate()
    )

    logging.info('Reading genes and initializeing liftover.')

    # Initialize LiftOver and gene objects:
    gene_index = spark.read.parquet(gene_index_file)
    lift = LiftOverSpark(chain_file)

    # Initialze the parser:
    logging.info('Starting Javierre data processing.')
    javierre = parser_javierre(javierre_data_file, gene_index, lift)

    # run QC:
    logging.info('Running QC on the intervals.')
    javierre.qc_intervals()

    # Save data:
    logging.info(f'Saving data to {output_file}.')
    javierre.save_parquet(output_file)


if __name__ == '__main__':
    # Parse arguments:

    # Initialize logger:

    # Initialize SparkSession:

    # Call main:
    main()
