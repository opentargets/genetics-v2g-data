from datetime import date
import os
import pandas as pd

import hydra

import pyspark
import pyspark.sql.types as T
import pyspark.sql.functions as F
from pyspark.sql import dataframe, SparkSession

from modules.andersson2014_parser import parse_anderson
from modules.Liftover import LiftOverSpark

# Get real path for the file:
script_folder = os.path.dirname(os.path.realpath(__file__))

@hydra.main(config_path=f'{script_folder}/configs', config_name='config')
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
    lift = LiftOverSpark(chain_file, max_difference)

    # Parsing the anderson file: <- Once more parsers are added a more elegant solution will be added.
    anderson_df = parse_anderson(cfg.intervals.anderson_file, gene_index, lift, proximity_limit).get_anderson_intervals()

    # Further parsers will come here...

    # Saving data:
    version = date.today().strftime("%y%m%d")
    anderson_df.write.mode('overwrite').parquet(cfg.intervals.output + f'/interval_{version}')


if __name__ == '__main__':

    # Calling parsers:
    main()
