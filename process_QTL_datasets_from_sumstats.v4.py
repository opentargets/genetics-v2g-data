#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
'''
This is the version of the script for version 4 (post- June 2019),
for version 3 (June 2019) I will hack the data so that the main
pipeline doesn't require any changes
'''



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
    
    # Repartition
    df = (
        df.repartitionByRange('chrom', 'pos')
        .orderBy('chrom', 'pos')
    )

    # Save
    (
        df
        .write
        .partitionBy('type', 'study_id', 'bio_feature')
        .json(
            outf,
            mode='overwrite'
        )
    )
    
    return 0

if __name__ == '__main__':

    main()
