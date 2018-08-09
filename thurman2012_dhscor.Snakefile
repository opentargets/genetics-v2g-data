#!/usr/bin/env snakemake
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import pandas as pd
from pprint import pprint

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']
UPLOAD = False

# Initiate remote handles
GS = GSRemoteProvider()

targets = []

### Make targets for Thurman 2012 DHS-promoter correlation dataset

# Processed
target = '{bucket}/{gs_dir}/{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.{proc}.tsv.gz'.format(
    bucket=config['gs_bucket'],
    gs_dir=config['gs_dir'],
    data_type='interval',
    exp_type='dhscor',
    source='thurman2012',
    cell_type='unspecified',
    proc='processed',
    chrom='1-23')
if UPLOAD:
    targets.append(GS.remote(target))

# Split files
for i in range(config['interval_split']):
    target = '{bucket}/{gs_dir}/{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.{proc}.split{i:03d}.tsv.gz'.format(
        bucket=config['gs_bucket'],
        gs_dir=config['gs_dir'],
        data_type='interval',
        exp_type='dhscor',
        source='thurman2012',
        cell_type='unspecified',
        proc='processed',
        chrom='1-23',
        i=i)
    if UPLOAD:
        targets.append(GS.remote(target))

# Raw
target = '{bucket}/{gs_dir}/{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.{proc}.bed.gz'.format(
    bucket=config['gs_bucket'],
    gs_dir=config['gs_dir'],
    data_type='interval',
    exp_type='dhscor',
    source='thurman2012',
    cell_type='unspecified',
    proc='raw',
    chrom='1-23')
if UPLOAD:
    targets.append(GS.remote(target))

# "all" must be first rule that is encounterd
rule all:
    ''' Master rule to trigger all targets (defined above) '''
    input:
        targets

# Import workflows
include: 'scripts/ensembl_grch37.Snakefile'

rule thurman2012_download:
    ''' Retrieves DHS correlation data from the Ensembl ENCODE FTP
    '''
    input:
        HTTPRemoteProvider().remote('http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/dhs_gene_connectivity/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz')
    output:
        tmpdir + '/interval/dhscor/thurman2012/unspecified/dhs_correlations.bed.gz'
    shell:
        'cp {input} {output}'

rule thurman2012_to_final:
    ''' Produces finalised output of fantom5 intervals
    '''
    input:
        bed = tmpdir + '/interval/dhscor/thurman2012/unspecified/dhs_correlations.bed.gz',
        gtf = tmpdir + '/Homo_sapiens.GRCh37.87.gtf.gz'
    output:
        tmpdir + '/interval/dhscor/thurman2012/unspecified/1-23.processed.tsv.gz'
    shell:
        'python scripts/thurman2012_to_final.py '
        '--inf {input.bed} '
        '--outf {output} '
        '--gtf {input.gtf} '
        '--cell_name Unspecified'

rule thurman2012_final_to_gcs:
    ''' Copy final to google cloud storage
    '''
    input:
        tmpdir + '/interval/dhscor/thurman2012/unspecified/1-23.processed.tsv.gz'
    output:
        GSRemoteProvider().remote(
            '{bucket}/{gs_dir}/interval/dhscor/thurman2012/unspecified/1-23.processed.tsv.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir']))
    shell:
        'cp {input} {output}'

rule thurman2012_raw_to_gcs:
    ''' Copy raw to google cloud storage
    '''
    input:
        tmpdir + '/interval/dhscor/thurman2012/unspecified/dhs_correlations.bed.gz'
    output:
        GSRemoteProvider().remote(
            '{bucket}/{gs_dir}/interval/dhscor/thurman2012/unspecified/1-23.raw.bed.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir']))
    shell:
        'cp {input} {output}'
