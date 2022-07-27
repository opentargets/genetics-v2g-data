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


class parse_thurman:
    """
    Parser Thurman 2012 dataset

    :param Thurman_parquet: path to the parquet file containing the Thurman 2016 data
    :param gene_index: Pyspark dataframe containing the gene index
    :param lift: LiftOverSpark object

    **Summary of the logic:**

    -
    """

    # Constants:
    DATASET_NAME = 'thurman2012'
    DATA_TYPE = 'interval'
    EXPERIMENT_TYPE = 'dhscor'
    PMID = '22955617'
    BIO_FEATURE = 'aggregate'

    def __init__(self,
                 thurman_datafile: str,
                 gene_index: dataframe,
                 lift: LiftOverSpark) -> None:

        logging.info('Parsing Thurman 2016 data...')

        thurman_schema = T.StructType([
            T.StructField('gene_chr', T.StringType(), False),
            T.StructField('gene_start', T.IntegerType(), False),
            T.StructField('gene_end', T.IntegerType(), False),
            T.StructField('gene_name', T.StringType(), False),
            T.StructField('chrom', T.StringType(), False),
            T.StructField('start', T.IntegerType(), False),
            T.StructField('end', T.IntegerType(), False),
            T.StructField('score', T.FloatType(), False)
        ])

        # Process Thurman data in a single step:
        self.Thurman_intervals = (
            SparkSession.getActiveSession()

            # Read table according to the schema, then do some modifications:
            .read.csv(thurman_datafile, sep='\t', header=False, schema=thurman_schema)
            .select(
                F.regexp_replace(F.col('chrom'), 'chr', '').alias('chrom'),
                'start', 'end', 'gene_name', 'score'
            )

            # Lift over to the GRCh38 build:
            .transform(lambda df: lift.convert_intervals(df, 'chrom', 'start', 'end'))

            # Map gene names to gene IDs:
            .join(gene_index.select('gene_id', 'gene_name'), how='inner', on='gene_name')

            # Select relevant columns and add constant columns:
            .select(
                'chrom',
                F.col('mapped_start').alias('start'),
                F.col('mapped_end').alias('end'),
                'gene_id', 'score',
                F.lit(self.DATASET_NAME).alias('dataset_name'),
                F.lit(self.DATA_TYPE).alias('data_type'),
                F.lit(self.EXPERIMENT_TYPE).alias('experiment_type'),
                F.lit(self.PMID).alias('pmid'),
                F.lit(self.BIO_FEATURE).alias('bio_feature'),
                F.lit(None).cast(T.StringType()).alias('cell_type'),
                F.lit(None).cast(T.StringType()).alias('tissue'),
            )
            .distinct()
            .persist()
        )

    def get_intervals(self) -> dataframe:
        return self.Thurman_intervals

    def qc_intervals(self) -> None:
        """
        Perform QC on the anderson intervals.
        """

        # Get numbers:
        logging.info(f'Size of Thurman data: {self.Thurman_intervals.count()}')
        logging.info(f'Number of unique intervals: {self.Thurman_intervals.select("start", "end").distinct().count()}')
        logging.info(f'Number genes in the Thurman dataset: {self.Thurman_intervals.select("gene_id").distinct().count()}')

    def save_parquet(self, output_file: str) -> None:
        self.Thurman_intervals.write.mode('overwrite').parquet(output_file)


def main(thurman_data_file: str, gene_index_file: str, chain_file: str, output_file: str) -> None:

    spark_conf = (
        SparkConf()
        .set('spark.driver.memory', '10g')
        .set('spark.executor.memory', '10g')
        .set('spark.driver.maxResultSize', '0')
        .set('spark.debug.maxToStringFields', '2000')
        .set('spark.sql.execution.arrow.maxRecordsPerBatch', '500000')
    )
    spark = (
        pyspark.sql.SparkSession.builder.config(conf=spark_conf)
        .master('local[*]')
        .config("spark.driver.bindAddress", "127.0.0.1")
        .getOrCreate()
    )

    logging.info('Reading genes and initializeing liftover.')

    # Initialize LiftOver and gene objects:
    gene_index = spark.read.parquet(gene_index_file)
    lift = LiftOverSpark(chain_file)

    # Initialze the parser:
    logging.info('Starting Thurman data processing.')
    thurman = parse_thurman(thurman_data_file, gene_index, lift)

    # run QC:
    logging.info('Running QC on the intervals.')
    thurman.qc_intervals()

    # Save data:
    logging.info(f'Saving data to {output_file}.')
    thurman.save_parquet(output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Wrapper for the the Thurman interval data parser.')
    parser.add_argument('--thurman_file', type=str, help='Path to the tsv dataset.')
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
    logging.info(f'Thurman file: {args.thurman_file}')
    logging.info(f'Gene index file: {args.gene_index}')
    logging.info(f'Chain file: {args.chain_file}')
    logging.info(f'Output file: {args.output_file}')

    main(
        args.thurman_file,
        args.gene_index,
        args.chain_file,
        args.output_file
    )
