from datetime import date
from functools import reduce
import os
import pandas as pd

import hydra

import pyspark
import pyspark.sql.types as T
import pyspark.sql.functions as F
from pyspark.sql import dataframe, SparkSession

from modules.andersson2014_parser import parse_anderson
from modules.javierre2016_parser import parse_javierre
from modules.jung2019_parser import parse_jung
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

    # Open and process gene file:
    gene_index = spark.read.parquet(cfg.intervals.gene_index).persist()

    # Initialize liftover object:
    lift = LiftOverSpark(chain_file, max_difference)

    # Parsing datasets:
    datasets = [

        # Parsing Andersson data:
        parse_anderson(cfg.intervals.anderson_file, gene_index, lift).get_intervals(),

        # Parsing Javierre data:
        parse_javierre(cfg.intervals.javierre_dataset, gene_index, lift).get_intervals(),

        # Parsing jung data:
        parse_jung(cfg.intervals.jung_dataset, gene_index, lift).get_intervals(),
    ]

    # Further parsers will come here...
    df = reduce(lambda x, y: x.union(y), datasets)

    # Saving data:
    version = date.today().strftime("%y%m%d")
    df.write.mode('overwrite').parquet(cfg.intervals.output + f'/interval_{version}')


if __name__ == '__main__':

    # Calling parsers:
    main()
