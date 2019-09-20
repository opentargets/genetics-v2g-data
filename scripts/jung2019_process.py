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

List of tissues:
['Trophoblast',                                                                 
 'Left Ventricle',
 'Adrenal Gland',
 'Human ES cells',
 'Neural Progenitor Cell',
 'Small Bowel',
 'Right Ventricle',
 'Mesendoderm',
 'Mesenchymal Stem Cell',
 'Spleen',
 'Sigmoid Colon',
 'Fibroblast',
 'Hippocampus',
 'Esophagus',
 'Dorsolateral Prefrontal Cortex',
 'Fat',
 'Bladder',
 'Thymus',
 'Liver',
 'Psoas',
 'Lymphoblastoid cell lines',
 'Pancreas',
 'Aorta',
 'Lung',
 'Gastric',
 'Right Atrium',
 'Ovary']

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

    #
    # Load --------------------------------------------------------------------
    #

    # Load Jung data
    import_schema = (
        StructType()
        .add('gene_names', StringType(), nullable=False)
        .add('interval', StringType(), nullable=False)
        .add('biofeature_string', StringType(), nullable=False)
    )
    data = spark.read.csv(
        args.indata,
        schema=import_schema,
        header=True,
        sep=','
    )

    # Load cell map
    cellmap = (
        spark.read.json(args.cellmap)
        .select(
            'biofeature_string', 'biofeature_code'
        )
        .distinct()
    )
    
    # Load gene data
    genes = (
        spark.read.json(args.genes)
    )

    #
    # Create gene symbol + synonym to gene id map -----------------------------
    #

    # Make official symbol map
    gene_official_map = (
        genes
        .select(
            'gene_name',
            'gene_id'
        )
    )

    # Make synonym map
    gene_synonym_map = (
        genes
        .select(
            explode(col('gene_synonyms')).alias('gene_name'),
            'gene_id'
        )
    )
    
    # Remove rows from the synonym table if they have an official match
    gene_synonym_map = (
        gene_synonym_map.join(
            gene_official_map,
            on='gene_name',
            how='left_anti'
        )
    )

    # Combine official and synonym tables
    gene_map = gene_official_map.unionByName(gene_synonym_map)
    
    #
    # Process data ------------------------------------------------------------
    #

    # Split, explode and map gene names to gene_id
    data = (
        data
        # Split and explode
        .withColumn('gene_name', explode(split(col('gene_names'), ';')))
        .drop('gene_names')
        # Add gene ids
        .join(
            gene_map,
            on='gene_name',
            how='left'
        )
        # Drop gene name
        .drop('gene_name')
    )

    # Split interval into chrom, start, end
    data = (
        data
        .withColumn('interval', split(col('interval'), '\.'))
        .withColumn('chrom_b37', col('interval').getItem(0))
        .withColumn('start_b37', col('interval').getItem(1).cast('long'))
        .withColumn('end_b37', col('interval').getItem(2).cast('long'))
        .drop('interval')
        # Strip "chr" from chomosomes
        .withColumn('chrom_b37',
            regexp_replace(col('chrom_b37'), 'chr', '')
        )
    )

    # Map tissues
    data = data.join(
        broadcast(cellmap),
        on='biofeature_string',
        how='left'
    ).cache()

    # Make sure there are no missing tissue codes
    assert data.filter(col('biofeature_code').isNull()).count() == 0

    #
    # Write -------------------------------------------------------------------
    #

    # Write to tsv using pandas (so that its compatible with the liftover script)
    (
        data
        .withColumn('score', lit(1))
        .select(
            col('chrom_b37').alias('chrom'),
            col('start_b37').alias('start'),
            col('end_b37').alias('end'),
            'gene_id',
            'score',
            col('biofeature_code').alias('bio_feature')
        )
        .toPandas()
        .to_csv(
            args.outdata,
            sep='\t',
            index=None,
            compression='gzip'
        )
    )

    return 0

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--indata', metavar="<file>", help=('Input file'), type=str, nargs='+', required=True)
    parser.add_argument('--genes', metavar="<file>", help=('Ensembl gene dictionary'), type=str, nargs='+', required=True)
    parser.add_argument('--cellmap', metavar="<file>", help=('Json to map tissue names to standard name'), type=str, nargs='+', required=True)
    parser.add_argument('--outdata', metavar="<file>", help=("Output file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
