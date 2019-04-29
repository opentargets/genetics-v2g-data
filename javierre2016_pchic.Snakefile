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

# Make targets for Javierre 2016 PCHiC

# Get list of cell line names
cell_types,  = GSRemoteProvider().glob_wildcards('gs://genetics-portal-input/v2g_input/javierre2016/{samples}.merged_samples_12Apr2015_full.txt.gz')

# # DEBUG
# print("WARNING! Only running 1 cell type")
# cell_types = cell_types[:1]

# Create output
target = '{out_dir}/{data_type}/{exp_type}/{source}/{version}/data.parquet'.format(
    out_dir=config['out_dir'],
    data_type='interval',
    exp_type='pchic',
    source='javierre2016',
    version=version)
targets.append(target)

# "all" must be first rule that is encounterd
rule all:
    ''' Master rule to trigger all targets (defined above) '''
    input:
        targets

# Import workflows
include: 'scripts/ensembl_grch37.Snakefile'

rule javierre2016_download:
    ''' Retrieves Javierre 2016 PCHiC file from the ebi FTP
    '''
    input:
        GSRemoteProvider().remote('gs://genetics-portal-input/v2g_input/javierre2016/{cell}.merged_samples_12Apr2015_full.txt.gz',
                                   keep_local=False, immediate_close=True)
    output:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/raw.gz'
    shell:
        'cp {input} {output}'

rule javierre2016_to_bed:
    ''' Formats Javierre file to bed.
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/raw.gz'
    output:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/raw.bed.gz'
    shell:
        'python scripts/javierre2016_to_bed.py '
        '--inf {input} '
        '--outf {output}'

rule javierre2016_tss_intersect:
    ''' Finds the intersect between PCHiC capture regions and gene TSS
    '''
    input:
        pchic = tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/raw.bed.gz',
        tss = tmpdir + '/Homo_sapiens.GRCh37.87.tss.protein_coding.bed.gz'
    output:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/tss.bed.gz'
    shell:
        'bedtools intersect -wa -wb '
        '-a {input.pchic} '
        '-b {input.tss} | '
        'gzip -c > {output}'

rule javierre2016_format:
    ''' Outputs formatted tsvs for Javierre dataset
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/tss.bed.gz'
    output:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/processed_b37.tsv.gz'
    shell:
        'python scripts/javierre2016_format.py '
        '--inf {input} '
        '--outf {output} '
        '--cell_name {wildcards.cell}'

rule liftover_to_GRCh38:
    ''' Lifts over co-ordinates to build 38
    '''
    input:
        data = tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/processed_b37.tsv.gz',
        chainfile = config['chain_grch37_to_grch38']
    output:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/processed_b38.tsv.gz'
    params:
        maxdiff = config['max_len_dff']
    shell:
        'python scripts/liftover_interval.py '
        '--inf {input.data} '
        '--outf {output} '
        '--chainfile {input.chainfile} '
        '--maxdiff {params.maxdiff} '

rule javierre2016_to_parquet:
    ''' Uses spark to write parquet file
    '''
    input:
        data=[tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/processed_b38.tsv.gz'.format(
            version=version, cell=cell_type) for cell_type in cell_types],
        cell_map = GSRemoteProvider().remote(config['cell_map'], keep_local=False, immediate_close=True)
    output:
        directory(config['out_dir'] + '/interval/pchic/javierre2016/{version}/data.parquet')
    shell:
        'python scripts/javierre2016_to_parquet.py '
        '--inf {input.data} '
        '--outf {output} '
        '--cell_map {input.cell_map} '
