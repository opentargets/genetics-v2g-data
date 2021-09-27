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
import pandas as pd

def main():

    # Args
    in_path = 'gs://genetics-portal-dev-sumstats/filtered/pvalue_0.05/molecular_trait/210917'
    outf = 'gs://genetics-portal-dev-staging/v2g/qtl/{version}'.format(
        version=date.today().strftime("%y%m%d")
    )

    # Make spark session
    global spark
    spark = (
        pyspark.sql.SparkSession.builder
        .getOrCreate()
    )
    print('Spark version: ', spark.version)
    
    # Load datasets
    df = spark.read.parquet(in_path)

    # Filter based on bonferonni correction of number of tests per gene
    df = df.filter(col('pval') <= (0.05 / col('num_tests')))

    # Only keep runs where gene_id is not null
    df = df.filter(col('gene_id').isNotNull())
    
    # Repartition
    df = (
        df.repartitionByRange('chrom', 'pos')
        .sortWithinPartitions('chrom', 'pos')
    ).persist()

    # Rename/select columns to what backend expects
    df = (
        df.withColumn('source', col('type_id'))
          .withColumnRenamed('type_id', 'type')
          .withColumnRenamed('gene_id', 'ensembl_id')
          .withColumnRenamed('ref', 'other_allele')
          .withColumnRenamed('alt', 'effect_allele')
          .select('type', 'source', 'study_id', 'bio_feature',
                  'chrom', 'pos', 'other_allele', 'effect_allele',
                  'ensembl_id', 'beta', 'se', 'pval')
    )

    # Save
    (
        df
        .write
        .partitionBy('type', 'study_id', 'bio_feature')
        .parquet(
            outf,
            mode='overwrite'
        )
    )
    
    # Save a list of all (source, features)
    (
        df
        .select('type', 'study_id', 'bio_feature')
        .drop_duplicates()
        .toPandas()
        .to_json(outf + '.feature_list.json', orient='records', lines=True)
    )

    return 0

if __name__ == '__main__':

    main()
