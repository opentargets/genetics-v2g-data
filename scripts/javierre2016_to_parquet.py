#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Ed Mountjoy
'''

'''
# Set SPARK_HOME and PYTHONPATH to use 2.4.0
export PYSPARK_SUBMIT_ARGS="--driver-memory 8g pyspark-shell"
export SPARK_HOME=/Users/em21/software/spark-2.4.0-bin-hadoop2.7
export PYTHONPATH=$SPARK_HOME/python:$SPARK_HOME/python/lib/py4j-2.4.0-src.zip:$PYTHONPATH
'''

import sys
import argparse
import pyspark.sql
from pyspark.sql.types import *
from pyspark.sql.functions import *
from functools import partial
import json

def main():

    # Parse args
    args = parse_args()

    # Make spark session
    spark = pyspark.sql.SparkSession.builder.getOrCreate()
    sc = spark.sparkContext
    print('Spark version: ', spark.version)

    # Make schema for import
    import_schema = (
        StructType()
        .add('chrom', StringType(), nullable=False)
        .add('start', LongType(), False)
        .add('end', LongType(), False)
        .add('gene_id', StringType(), False)
        .add('score', DoubleType(), False)
        .add('bio_feature', StringType(), False)
    )

    # Load
    df = spark.read.csv(
        args.infs,
        schema=import_schema,
        header=True,
        sep='\t'
    )

    # Show schema
    df.printSchema()

    # Load cell map and broadcast
    cell_map_dict = sc.broadcast(load_cell_map(args.cell_map))

    # Map bio_feature string to code
    cell_map_func_udf = udf(
        partial(cell_map_func, d=cell_map_dict), StringType() )
    df = (
        df.withColumn('bio_feature_code', cell_map_func_udf(df.bio_feature))
          .drop('bio_feature')
          .withColumnRenamed('bio_feature_code', 'bio_feature') )

    # Repartition
    df = (
        df.repartitionByRange('chrom', 'start', 'end')
        .sortWithinPartitions('chrom', 'start', 'end')
    )
    print('Num partitions: ', df.rdd.getNumPartitions())

    # Write
    ( df.write.parquet(
          args.outf,
          mode='overwrite',
          compression='snappy')
    )

    return 0

def cell_map_func(x, d):
    return d.value[x]

def load_cell_map(inf):
    d = {}
    with open(inf, 'r') as in_h:
        for line in in_h:
            parts = json.loads(line)
            d[parts['biofeature_string']] = parts['biofeature_code']
    return d

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--infs', metavar="<file>", help=('Input file'), type=str, nargs='+', required=True)
    parser.add_argument('--outf', metavar="<file>", help=("Output file"), type=str, required=True)
    parser.add_argument('--cell_map', metavar="<file>", help=("Cell amp"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
