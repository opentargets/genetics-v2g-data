#!/usr/bin/env snakemake
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
import pandas as pd
from pprint import pprint
from datetime import date

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']
version = date.today().strftime("%y%m%d")

targets = []

### Make targets for Thurman 2012 DHS-promoter correlation dataset

# Processed
target = '{out_dir}/{data_type}/{exp_type}/{source}/{version}/data.parquet'.format(
    out_dir=config['out_dir'],
    data_type='interval',
    exp_type='dhscor',
    source='thurman2012',
    version=version,
    cell_type='aggregate')
targets.append(target)

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
        GSRemoteProvider().remote('gs://genetics-portal-input/v2g_input/thurman2012/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz',
                                    keep_local=False)
    output:
        tmpdir + '/interval/dhscor/thurman2012/unspecified/dhs_correlations.bed.gz'
    shell:
        'cp {input} {output}'

rule thurman2012_format:
    ''' Produces formatted output of fantom5 intervals
    '''
    input:
        bed = tmpdir + '/interval/dhscor/thurman2012/unspecified/dhs_correlations.bed.gz',
        gtf = tmpdir + '/Homo_sapiens.GRCh37.87.gtf.gz'
    output:
        tmpdir + '/interval/dhscor/thurman2012/{version}/data_b37.tsv.gz'
    shell:
        'python scripts/thurman2012_format.py '
        '--inf {input.bed} '
        '--outf {output} '
        '--gtf {input.gtf} '
        '--cell_name aggregate'

rule liftover_to_GRCh38:
    ''' Lifts over co-ordinates to build 38
    '''
    input:
        data = tmpdir + '/interval/dhscor/thurman2012/{version}/data_b37.tsv.gz',
        chainfile = config['chain_grch37_to_grch38']
    output:
        tmpdir + '/interval/dhscor/thurman2012/{version}/data_b38.tsv.gz'
    params:
        maxdiff = config['max_len_dff']
    shell:
        'python scripts/liftover_interval.py '
        '--inf {input.data} '
        '--outf {output} '
        '--chainfile {input.chainfile} '
        '--maxdiff {params.maxdiff} '

rule thurman2012_to_parquet:
    ''' Uses spark to write parquet file
    '''
    input:
        tmpdir + '/interval/dhscor/thurman2012/{version}/data_b38.tsv.gz'
    output:
        directory(config['out_dir'] + '/interval/dhscor/thurman2012/{version}/data.parquet')
    shell:
        'python scripts/thurman2012_to_parquet.py '
        '--inf {input} '
        '--outf {output}'
