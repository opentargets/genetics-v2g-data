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
    in_path = 'gs://genetics-portal-sumstats-b38/unfiltered/molecular_trait'
    study_list = [
        'ALASOO_2018.parquet',
        'Blueprint.parquet',
        'CEDAR.parquet',
        'FAIRFAX_2012.parquet',
        'FAIRFAX_2014.parquet',
        'GENCORD.parquet',
        'GEUVADIS.parquet',
        'GTEX_v7.parquet',
        'HIPSCI.parquet',
        'NARANBHAI_2015.parquet',
        'NEDELEC_2016.parquet',
        'QUACH_2016.parquet',
        'SCHWARTZENTRUBER_2018.parquet',
        'SUN2018.parquet',
        'TWINSUK.parquet',
        'VAN_DE_BUNT_2015.parquet',
        'eQTLGen.parquet'
    ]
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

    # Load list of datasets
    dfs = []
    for in_path in [os.path.join(in_path, study) for study in study_list]:
        df_temp = spark.read.parquet(in_path)
        dfs.append(df_temp)
    
    # Take union
    df = reduce(pyspark.sql.DataFrame.unionByName, dfs)

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
