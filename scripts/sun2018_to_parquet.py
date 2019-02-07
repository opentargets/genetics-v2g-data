#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Ed Mountjoy
'''

import sys
import argparse
import pyspark.sql
from pyspark.sql.types import *
from pyspark.sql.functions import *
from functools import partial

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
        .add('pos', LongType(), False)
        .add('ref', StringType(), False)
        .add('alt', StringType(), False)
        .add('beta', DoubleType(), False)
        .add('se', DoubleType(), False)
        .add('pval', DoubleType(), False)
    )

    # Load
    df = spark.read.csv(
        args.infs,
        schema=import_schema,
        header=True,
        sep='\t'
    )

    # Extract add bio_feaute
    get_gene_udf = udf(get_gene)
    df = df.withColumn('gene_id', get_gene_udf(input_file_name()))

    # Add bio_feature
    df = df.withColumn('bio_feaute', lit(args.cell_code))

    # Rearrange columns
    df = df.select('chrom', 'pos', 'ref', 'alt', 'gene_id', 'beta', 'se', 'pval', 'bio_feaute')

    # df.show(3)
    # df.printSchema()

    # Repartition
    df = df.repartitionByRange('chrom', 'pos', 'ref', 'alt')
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

def get_gene(filename):
    ''' Returns gene from filename
    '''
    return filename.split('/')[-1].split('.')[0]

def load_cell_map(inf):
    d = {}
    with open(inf, 'r') as in_h:
        in_h.readline()
        for line in in_h:
            original, label, code = line.rstrip().split(',')
            d[original] = code
    return d

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--infs', metavar="<file>", help=('Input file'), type=str, nargs='+', required=True)
    parser.add_argument('--outf', metavar="<file>", help=("Output file"), type=str, required=True)
    parser.add_argument('--cell_code', metavar="<file>", help=("Cell code"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
