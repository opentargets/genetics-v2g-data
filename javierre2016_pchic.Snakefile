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

# Make targets for Javierre 2016 PCHiC

Get list of cell line names
cell_types,  = FTPRemoteProvider().glob_wildcards('ftp.ebi.ac.uk/pub/contrib/pchic/'
    'CHiCAGO/{samples}.merged_samples_12Apr2015_full.txt.gz')

# Add target file for each cell line
for cell_type in list(cell_types):

    # Processed
    target = '{bucket}/{gs_dir}/{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.{proc}.tsv.gz'.format(
        bucket=config['gs_bucket'],
        gs_dir=config['gs_dir'],
        data_type='interval',
        exp_type='pchic',
        source='javierre2016',
        cell_type=cell_type,
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
            exp_type='pchic',
            source='javierre2016',
            cell_type=cell_type,
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
        exp_type='pchic',
        source='javierre2016',
        cell_type=cell_type,
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

rule javierre2016_download:
    ''' Retrieves Javierre 2016 PCHiC file from the ebi FTP
    '''
    input:
        FTPRemoteProvider().remote('ftp.ebi.ac.uk/pub/contrib/pchic/CHiCAGO/{cell}.merged_samples_12Apr2015_full.txt.gz')
    output:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.raw.gz'
    shell:
        'cp {input} {output}'

rule javierre2016_to_bed:
    ''' Formats Javierre file to bed.
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.raw.gz'
    output:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.raw.bed.gz'
    shell:
        'python scripts/javierre2016_to_bed.py '
        '--inf {input} '
        '--outf {output}'

rule javierre2016_tss_intersect:
    ''' Finds the intersect between PCHiC capture regions and gene TSS
    '''
    input:
        pchic = tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.raw.bed.gz',
        tss = tmpdir + '/Homo_sapiens.GRCh37.87.tss.bed.gz'
    output:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.tss.bed.gz'
    shell:
        'bedtools intersect -wa -wb '
        '-a {input.pchic} '
        '-b {input.tss} | '
        'gzip -c > {output}'

rule javierre2016_to_final:
    ''' Outputs finalised format for Javierre dataset
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.tss.bed.gz'
    output:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.final.tsv.gz'
    shell:
        'python scripts/javierre2016_to_final.py '
        '--inf {input} '
        '--outf {output} '
        '--cell_name {wildcards.cell}'

rule javierre2016_final_to_gcs:
    ''' Copy final to google cloud storage
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.final.tsv.gz'
    output:
        GSRemoteProvider().remote(
            '{bucket}/{gs_dir}/interval/pchic/javierre2016/{{cell}}/1-23.processed.tsv.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir']))
    shell:
        'cp {input} {output}'

rule javierre2016_raw_to_gcs:
    ''' Copy raw to google cloud storage
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.raw.gz'
    output:
        GSRemoteProvider().remote(
            '{bucket}/{gs_dir}/interval/pchic/javierre2016/{{cell}}/1-23.raw.bed.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir']))
    shell:
        'cp {input} {output}'
