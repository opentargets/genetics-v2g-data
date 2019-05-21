#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

'''
# Set SPARK_HOME and PYTHONPATH to use 2.4.0
export PYSPARK_SUBMIT_ARGS="--driver-memory 8g pyspark-shell"
export SPARK_HOME=/Users/em21/software/spark-2.4.0-bin-hadoop2.7
export PYTHONPATH=$SPARK_HOME/python:$SPARK_HOME/python/lib/py4j-2.4.0-src.zip:$PYTHONPATH
'''

import os
import sys
import pyspark.sql
from pyspark.sql.types import *
from pyspark.sql.functions import *
from functools import reduce
from datetime import date

def main():

    # Args
    in_path = 'gs://genetics-portal-sumstats-b38/filtered/pvalue_0.05/molecular_trait'
    outf = 'gs://genetics-portal-staging/v2g/qtl/{version}'.format(
        version=date.today().strftime("%y%m%d")
    )

    # # Args (test)
    # in_path = 'tmp/part-00000-1a32ece0-aba4-42f2-b9b6-963e49610fba-c000.json'
    # outf = 'tmp/qtl/{version}'.format(
    #     version=date.today().strftime("%y%m%d")
    # )

    # Make spark session
    global spark
    spark = (
        pyspark.sql.SparkSession.builder
        .getOrCreate()
    )
    print('Spark version: ', spark.version)

    # Load datasets
    df = spark.read.json(in_path)

    # Filter based on bonferonni correction of number of tests per gene
    df = df.filter(col('pval') <= (0.05 / col('num_tests')))

    # Only keep runs where gene_id is not null
    df = df.filter(col('gene_id').isNotNull())

    # Hack the columns to fit the current structure
    hack = (
        df
        .withColumn('source', col('type'))
        .withColumn('feature', concat_ws('-', col('study_id'), col('bio_feature')))
        .withColumnRenamed('gene_id', 'ensembl_id')
        .withColumnRenamed('ref', 'other_allele')
        .withColumnRenamed('alt', 'effect_allele')
        .select('type', 'source', 'feature', 'chrom',
                'pos', 'other_allele', 'effect_allele',
                'ensembl_id', 'beta', 'se', 'pval'
        )
    )

    # Repartition
    hack = (
        hack.repartitionByRange('chrom', 'pos')
        .orderBy('chrom', 'pos')
    ).persist()

    # Save data
    (
        hack
        .write
        .parquet(
            outf,
            mode='overwrite'
        )
    )

    # Save a list of all (source, features)
    (
        hack
        .select('source', 'feature')
        .drop_duplicates()
        .coalesce(1)
        .write
        .json(
            outf + '.feature_list.json',
            mode='overwrite'
        )
    )


    
    return 0

if __name__ == '__main__':

    main()
