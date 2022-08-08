import argparse
from functools import reduce
import logging

from psutil import virtual_memory
import pyspark
import pyspark.sql.types as T
import pyspark.sql.functions as F
from pyspark.sql import dataframe, SparkSession
from pyspark.conf import SparkConf

from modules.Liftover import LiftOverSpark


class parse_javierre:
    """
    Parser Javierre 2016 dataset

    :param javierre_parquet: path to the parquet file containing the Javierre 2016 data after processing it (see notebooks/Javierre_data_pre-process.ipynb)
    :param gene_index: Pyspark dataframe containing the gene index
    :param lift: LiftOverSpark object

    **Summary of the logic:**

    - Reading parquet file containing the pre-processed Javierre dataset.
    - Splitting name column into chromosome, start, end, and score.
    - Lifting over the intervals.
    - Mapping intervals to genes by overlapping regions.
    - For each gene/interval pair, keep only the highest scoring interval.
    - Filter gene/interval pairs by the distance between the TSS and the start of the interval.
    """

    # Constants:
    DATASET_NAME = 'javierre2016'
    DATA_TYPE = 'interval'
    EXPERIMENT_TYPE = 'pchic'
    PMID = '27863249'
    TWOSIDED_THRESHOLD = 2.45e6  # <-  this needs to phased out. Filter by percentile instead of absolute value.

    def __init__(self,
                 javierre_parquet: str,
                 gene_index: dataframe,
                 lift: LiftOverSpark) -> None:

        logging.info('Parsing Javierre 2016 data...')

        # Read gene index:
        genes = self.prepare_gene_index(gene_index)

        # Read Javierre data:
        javierre_raw = (
            SparkSession.getActiveSession().read.parquet(javierre_parquet)

            # Splitting name column into chromosome, start, end, and score:
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
                (F.col('name_score').isNotNull()) &  # Dropping rows without score
                (F.col('chrom') == F.col('name_chr')) &  # 6_473_184 -> 6_428_824
                F.col('name_chr').isin([f'{x}' for x in range(1, 23)] + ['X', 'Y', 'MT'])  # Keep only canonical chromosomes
            )
        )

        # Lifting over intervals:
        javierre_remapped = (
            javierre_raw
            # Lifting over to GRCh38 interval 1:
            .transform(lambda df: lift.convert_intervals(df, 'chrom', 'start', 'end'))
            .drop('start', 'end')
            .withColumnRenamed('mapped_chrom', 'chrom')
            .withColumnRenamed('mapped_start', 'start')
            .withColumnRenamed('mapped_end', 'end')

            # Lifting over interval 2 to GRCh38:
            .transform(lambda df: lift.convert_intervals(df, 'name_chr', 'name_start', 'name_end'))
            .drop('name_start', 'name_end')
            .withColumnRenamed('mapped_name_chr', 'name_chr')
            .withColumnRenamed('mapped_name_start', 'name_start')
            .withColumnRenamed('mapped_name_end', 'name_end')
            .persist()
        )

        # Once the intervals are lifted, extracting the unique intervals:
        unique_intervals_with_genes = (
            javierre_remapped
            .select(
                'chrom',
                F.col('start').cast(T.IntegerType()),
                F.col('end').cast(T.IntegerType())
            )
            .distinct()
            .join(genes, on=['chrom'], how='left')
            .filter(
                ((F.col('start') >= F.col('gene_start')) & (F.col('start') <= F.col('gene_end'))) |
                ((F.col('end') >= F.col('gene_start')) & (F.col('end') <= F.col('gene_end')))
            )
            .select('chrom', 'start', 'end', 'gene_id', 'TSS')
        )

        # Joining back the data:
        self.javierre_intervals = (
            javierre_remapped
            .join(unique_intervals_with_genes, on=['chrom', 'start', 'end'], how='left')
            .filter(
                # Drop rows where the TSS is far from the start of the region
                F.abs((F.col('start') + F.col('end')) / 2 - F.col('TSS')) <= self.TWOSIDED_THRESHOLD
            )

            # For each gene, keep only the highest scoring interval:
            .groupBy('name_chr', 'name_start', 'name_end', 'gene_id', 'bio_feature')
            .agg(F.max(F.col('name_score')).alias('score'))

            # Create the output:
            .select(
                F.col('name_chr').alias('chrom'),
                F.col('name_start').alias('start'),
                F.col('name_end').alias('end'),
                F.col('score'),
                F.col('gene_id').alias('gene_id'),
                F.col('bio_feature').alias('cell_type'),
                F.lit(self.DATASET_NAME).alias('dataset_name'),
                F.lit(self.DATA_TYPE).alias('data_type'),
                F.lit(self.EXPERIMENT_TYPE).alias('experiment_type'),
                F.lit(self.PMID).alias('pmid')
            )
            .persist()
        )

    def get_intervals(self) -> dataframe:
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
        '''Pre-processing the gene dataset
        - selecting and renaming relevant columns
        - remove 'chr' from chromosome column

        :param gene_index: Path to the gene parquet file
        :return: Spark Dataframe
        '''
        # Reading gene annotations:
        return (
            gene_index
            .select(
                F.regexp_replace(F.col('chr'), 'chr', '').alias('chrom'),
                F.col('start').cast(T.IntegerType()).alias('gene_start'),
                F.col('end').cast(T.IntegerType()).alias('gene_end'),
                'gene_id', 'TSS'
            )
            .persist()
        )


def main(javierre_data_file: str, gene_index_file: str, chain_file: str, output_file: str) -> None:

    spark_conf = (
        SparkConf()
        .set('spark.driver.memory', '10g')
        .set('spark.executor.memory', '10g')
        .set('spark.driver.maxResultSize', '0')
        .set('spark.debug.maxToStringFields', '2000')
        .set('spark.sql.execution.arrow.maxRecordsPerBatch', '500000')
        .set('spark.driver.bindAddress', '127.0.0.1')
    )
    spark = (
        pyspark.sql.SparkSession.builder.config(conf=spark_conf)
        .master('local[*]')
        .getOrCreate()
    )

    logging.info('Reading genes and initializeing liftover.')

    # Initialize LiftOver and gene objects:
    gene_index = spark.read.parquet(gene_index_file)
    lift = LiftOverSpark(chain_file)

    # Initialze the parser:
    logging.info('Starting Javierre data processing.')
    javierre = parse_javierre(javierre_data_file, gene_index, lift)

    # run QC:
    logging.info('Running QC on the intervals.')
    javierre.qc_intervals()

    # Save data:
    logging.info(f'Saving data to {output_file}.')
    javierre.save_parquet(output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Wrapper for the the Javierre interval data parser.')
    parser.add_argument('--javierre_file', type=str, help='Path to the pre-processed parquet dataset (.parquet).')
    parser.add_argument('--gene_index', type=str, help='Path to the gene index file (.parquet)')
    parser.add_argument('--chain_file', type=str, help='Path to the chain file (.chain)')
    parser.add_argument('--output_file', type=str, help='Path to the output file (.parquet)')
    args = parser.parse_args()

    # Initialize logging:
    logging.basicConfig(
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    # Just print out some of the arguments:
    logging.info(f'Javierre file: {args.javierre_file}')
    logging.info(f'Gene index file: {args.gene_index}')
    logging.info(f'Chain file: {args.chain_file}')
    logging.info(f'Output file: {args.output_file}')

    main(
        args.javierre_file,
        args.gene_index,
        args.chain_file,
        args.output_file
    )
