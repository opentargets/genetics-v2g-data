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

### Make targets for Andersson 2014 Enhancer-TSS association dataset (Fantom5)

# Processed
target = '{out_dir}/{data_type}/{exp_type}/{source}/{version}/data.parquet'.format(
    out_dir=config['out_dir'],
    data_type='interval',
    exp_type='fantom5',
    source='andersson2014',
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

rule andersson2014_download:
    ''' Retrieves Andersson et al 2014 Fantom5
    '''
    input:
        GSRemoteProvider().remote('gs://genetics-portal-input/v2g_input/andersson2014/enhancer_tss_associations.bed',
                                    keep_local=False)
    output:
        tmpdir + '/interval/fantom5/andersson2014/aggregate/enhancer_tss_associations.bed'
    shell:
        'cp {input} {output}'

rule andersson2014_format:
    ''' Produces format output of fantom5 intervals
    '''
    input:
        bed = tmpdir + '/interval/fantom5/andersson2014/aggregate/enhancer_tss_associations.bed',
        gtf = GSRemoteProvider().remote('gs://genetics-portal-data/lut/gene_dictionary.json',
                                        keep_local=False)
    output:
        tmpdir + '/interval/fantom5/andersson2014/{version}/data_b37.tsv.gz'
    shell:
        'python scripts/andersson2014_format.py '
        '--inf {input.bed} '
        '--outf {output} '
        '--gene_info {input.gtf} '
        '--cell_name aggregate'

rule liftover_to_GRCh38:
    ''' Lifts over co-ordinates to build 38
    '''
    input:
        data = tmpdir + '/interval/fantom5/andersson2014/{version}/data_b37.tsv.gz',
        chainfile = config['chain_grch37_to_grch38']
    output:
        tmpdir + '/interval/fantom5/andersson2014/{version}/data_b38.tsv.gz'
    params:
        maxdiff = config['max_len_dff']
    shell:
        'python scripts/liftover_interval.py '
        '--inf {input.data} '
        '--outf {output} '
        '--chainfile {input.chainfile} '
        '--maxdiff {params.maxdiff} '

rule andersson2014_to_parquet:
    ''' Uses spark to write parquet file
    '''
    input:
        tmpdir + '/interval/fantom5/andersson2014/{version}/data_b38.tsv.gz'
    output:
        directory(config['out_dir'] + '/interval/fantom5/andersson2014/{version}/data.parquet')
    shell:
        'python scripts/andersson2014_to_parquet.py '
        '--inf {input} '
        '--outf {output}'
