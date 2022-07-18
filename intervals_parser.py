from datetime import date
import pandas as pd

import hydra

import pyspark.sql.types as T
import pyspark.sql.functions as F
from pyspark.sql import dataframe, SparkSession

from modules import andersson2014_parser, Liftover


@hydra.main(config_path='../configs', config_name='config')
def main(cfg):

    spark = (
        pyspark.sql.SparkSession
        .builder
        .master("local[*]")
        .getOrCreate()
    )

    chain_file = cfg.intervals.liftover_chain_file
    max_difference = cfg.intervals.max_length_difference
    proximity_limit = 2 * float(cfg.intervals.proximity_limit)

    # Open and process gene file:
    gene_index = spark.read.parquet(cfg.intervals.gene_index).persist()

    # Initialize liftover object:
    lift = Liftover.LiftOverSpark(chain_file, max_difference)

    # Parsing the anderson file:
    anderson_df = andersson2014_parser.parse_anderson(cfg.intervals.anderson, gene_index, lift, proximity_limit)

    # Further parsers will come here...

    # Saving data:
    version = date.today().strftime("%y%m%d")
    anderson_df.write.mode('overwrite').parquet(cfg.intervals.output + f'/interval_{version}')


if __name__ == '__main__':

    # Calling parsers:
    main()
