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

### Make targets for Andersson 2014 Enhancer-TSS association dataset (Fantom5)

# Processed
target = '{bucket}/{gs_dir}/{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.{proc}.tsv.gz'.format(
    bucket=config['gs_bucket'],
    gs_dir=config['gs_dir'],
    data_type='interval',
    exp_type='fantom5',
    source='andersson2014',
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
        exp_type='fantom5',
        source='andersson2014',
        cell_type='unspecified',
        proc='processed',
        chrom='1-23',
        i=i)
    if UPLOAD:
        targets.append(GS.remote(target))

# Raw
target = '{bucket}/{gs_dir}/{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.{proc}.bed'.format(
    bucket=config['gs_bucket'],
    gs_dir=config['gs_dir'],
    data_type='interval',
    exp_type='fantom5',
    source='andersson2014',
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

rule andersson2014_download:
    ''' Retrieves Andersson et al 2014 Fantom5
    '''
    input:
        HTTPRemoteProvider().remote('http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed')
    output:
        tmpdir + '/interval/fantom5/andersson2014/unspecified/enhancer_tss_associations.bed'
    shell:
        'cp {input} {output}'

rule andersson2014_to_final:
    ''' Produces finalised output of fantom5 intervals
    '''
    input:
        bed = tmpdir + '/interval/fantom5/andersson2014/unspecified/enhancer_tss_associations.bed',
        gtf = tmpdir + '/Homo_sapiens.GRCh37.87.gtf.gz'
    output:
        tmpdir + '/interval/fantom5/andersson2014/unspecified/1-23.processed.tsv.gz'
    shell:
        'python scripts/andersson2014_to_final.py '
        '--inf {input.bed} '
        '--outf {output} '
        '--gtf {input.gtf} '
        '--cell_name Unspecified'

rule andersson2014_final_to_gcs:
    ''' Copy to google cloud storage
    '''
    input:
        tmpdir + '/interval/fantom5/andersson2014/unspecified/1-23.processed.tsv.gz'
    output:
        GSRemoteProvider().remote(
            '{bucket}/{gs_dir}/interval/fantom5/andersson2014/unspecified/1-23.processed.tsv.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir']))
    shell:
        'cp {input} {output}'

rule andersson2014_raw_to_gcs:
    ''' Copy to google cloud storage
    '''
    input:
        tmpdir + '/interval/fantom5/andersson2014/unspecified/enhancer_tss_associations.bed'
    output:
        GSRemoteProvider().remote(
            '{bucket}/{gs_dir}/interval/fantom5/andersson2014/unspecified/1-23.raw.bed'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir']))
    shell:
        'cp {input} {output}'
