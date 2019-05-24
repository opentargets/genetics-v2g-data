#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Ed Mountjoy
'''

import sys
import argparse
import pyspark.sql
from pyspark.sql.types import *

def main():

    # Parse args
    args = parse_args()

    # Make spark session
    spark = pyspark.sql.SparkSession.builder.getOrCreate()
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
        args.inf,
        schema=import_schema,
        header=True,
        sep='\t'
    )

    # Show schema
    df.printSchema()

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

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('Input file'), type=str, required=True)
    parser.add_argument('--outf', metavar="<file>", help=("Output file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
