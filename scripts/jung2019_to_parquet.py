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

    # Load data
    data = spark.read.csv(
        args.inf,
        schema=import_schema,
        header=True,
        sep='\t'
    )

    # Load gene positions
    genes = (
        spark.read.json(args.genes)
        .select(
            'gene_id',
            'gene_chrom',
            'gene_start',
            'gene_end'
        )
        .withColumn('gene_mid',
            (col('gene_start') + col('gene_end')) / 2
        )
        .drop('gene_start', 'gene_end')
    )

    #
    # Remove inter-chromosomal and gene-interval distances > 2 sd -------------
    #

    # Add gene position information
    data = data.join(
        genes,
        on='gene_id',
        how='left'
    )

    # Remove inter-chromosomal interactions and calculate distance
    data = (
        data
        .filter(col('chrom') == col('gene_chrom'))
        .withColumn('distance',
            abs(((col('start') + col('end')) / 2) - col('gene_mid'))
        )
        .drop('gene_chrom', 'gene_mid')
    ).cache()

    # Calculate distance standard deviation
    dist_sd = (
        data
        .select(stddev(col('distance')))
        .rdd.map(lambda x: x[0])
        .collect()[0]
    )
    
    # Filter to only keep intervals < 2 * sd
    data = (
        data
        .filter(col('distance') <= 2 * dist_sd)
        .drop('distance')
    )

    #
    # Write -------------------------------------------------------------------
    #

    # Write
    ( 
        data
        .repartitionByRange('chrom', 'start', 'end')
        .sortWithinPartitions('chrom', 'start', 'end')
        .write.parquet(
            args.outf,
            mode='overwrite',
            compression='snappy'
        )
    )

    return 0

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('Input file'), type=str, required=True)
    parser.add_argument('--genes', metavar="<file>", help=('Gene info file'), type=str, required=True)
    parser.add_argument('--outf', metavar="<file>", help=("Output file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
