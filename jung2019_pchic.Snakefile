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

# Make targets for Jung 2019 PCHiC

# Create output
target = '{out_dir}/{data_type}/{exp_type}/{source}/{version}/data.parquet'.format(
    out_dir=config['out_dir'],
    data_type='interval',
    exp_type='pchic',
    source='jung2019',
    version=version)
targets.append(target)

# "all" must be first rule that is encounterd
rule all:
    ''' Master rule to trigger all targets (defined above) '''
    input:
        targets

rule javierre2016_process:
    ''' Processes Jung 2019 PCHiC to a GRCh37 interval file
    '''
    input:
        data = GSRemoteProvider().remote(
            'genetics-portal-raw/pchic_jung2019/jung2019_pchic_tableS3.csv',
            keep_local=True,
            immediate_close=True
        ),
        genes = GSRemoteProvider().remote(
            config['genes'],
            keep_local=True,
            immediate_close=True
        ),
        tissues = GSRemoteProvider().remote(
            config['cell_map'],
            keep_local=True,
            immediate_close=True
        )
    output:
        tmpdir + '/interval/pchic/jung2019/{version}/processed_b37.tsv.gz'
    shell:
        'python scripts/jung2019_process.py '
        '--indata {input.data} '
        '--genes {input.genes} '
        '--cellmap {input.tissues} '
        '--outdata {output}'

rule liftover_to_GRCh38:
    ''' Lifts over co-ordinates to build 38
    '''
    input:
        data = tmpdir + '/interval/pchic/jung2019/{version}/processed_b37.tsv.gz',
        chainfile = config['chain_grch37_to_grch38']
    output:
        tmpdir + '/interval/pchic/jung2019/{version}/processed_b38.tsv.gz'
    params:
        maxdiff = config['max_len_dff']
    shell:
        'python scripts/liftover_interval.py '
        '--inf {input.data} '
        '--outf {output} '
        '--chainfile {input.chainfile} '
        '--maxdiff {params.maxdiff} '

rule jung2019_to_parquet:
    ''' Uses spark to write parquet file
    '''
    input:
        data = tmpdir + '/interval/pchic/jung2019/{version}/processed_b38.tsv.gz',
        genes = GSRemoteProvider().remote(
            config['genes'],
            keep_local=True,
            immediate_close=True
        )
    output:
        directory(config['out_dir'] + '/interval/pchic/jung2019/{version}/data.parquet')
    shell:
        'python scripts/jung2019_to_parquet.py '
        '--inf {input.data} '
        '--genes {input.genes} '
        '--outf {output}'
