#!/usr/bin/env snakemake
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
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
cell_types,  = FTPRemoteProvider().glob_wildcards('ftp.ebi.ac.uk/pub/contrib/pchic/'
    'CHiCAGO/{samples}.merged_samples_12Apr2015_full.txt.gz')

# Add target file for each cell line
for cell_type in list(cell_types):

    # Processed
    target = '{out_dir}/{data_type}/{exp_type}/{source}/{version}/{cell_type}/{chrom}.{proc}.tsv.gz'.format(
        out_dir=config['out_dir'],
        data_type='interval',
        exp_type='pchic',
        source='javierre2016',
        version=version,
        cell_type=cell_type,
        proc='processed',
        chrom='1-23')
    targets.append(target)

    # Split files
    for i in range(config['interval_split']):
        target = '{out_dir}/{data_type}/{exp_type}/{source}/{version}/{cell_type}/{chrom}.{proc}.split{i:03d}.tsv.gz'.format(
            out_dir=config['out_dir'],
            data_type='interval',
            exp_type='pchic',
            source='javierre2016',
            version=version,
            cell_type=cell_type,
            proc='processed',
            chrom='1-23',
            i=i)
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
        FTPRemoteProvider().remote('ftp.ebi.ac.uk/pub/contrib/pchic/CHiCAGO/{cell}.merged_samples_12Apr2015_full.txt.gz',
                                   keep_local=False, immediate_close=True)
    output:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/1-23.raw.gz'
    shell:
        'cp {input} {output}'

rule javierre2016_to_bed:
    ''' Formats Javierre file to bed.
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/1-23.raw.gz'
    output:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/1-23.raw.bed.gz'
    shell:
        'python scripts/javierre2016_to_bed.py '
        '--inf {input} '
        '--outf {output}'

rule javierre2016_tss_intersect:
    ''' Finds the intersect between PCHiC capture regions and gene TSS
    '''
    input:
        pchic = tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/1-23.raw.bed.gz',
        tss = tmpdir + '/Homo_sapiens.GRCh37.87.tss.protein_coding.bed.gz'
    output:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/1-23.tss.bed.gz'
    shell:
        'bedtools intersect -wa -wb '
        '-a {input.pchic} '
        '-b {input.tss} | '
        'gzip -c > {output}'

rule javierre2016_to_final:
    ''' Outputs finalised format for Javierre dataset
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/1-23.tss.bed.gz'
    output:
        config['out_dir'] + '/interval/pchic/javierre2016/{version}/{cell}/1-23.processed.tsv.gz'
    shell:
        'python scripts/javierre2016_to_final.py '
        '--inf {input} '
        '--outf {output} '
        '--cell_name {wildcards.cell}'

rule split_unzip_final:
    ''' Unzip and remove header in preparation for splitting
    '''
    input:
        config['out_dir'] + '/interval/pchic/javierre2016/{version}/{cell}/1-23.processed.tsv.gz'
    output:
        temp(tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/1-23.processed.tsv')
    shell:
        'zcat < {input} | tail -n +2 > {output}'

rule split:
    ''' Splits file into many parts
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/1-23.processed.tsv'
    output:
        temp(expand(tmpdir + '/interval/pchic/javierre2016/{{version}}/{{cell}}/1-23.processed.split{i:03d}.tsv',
               i=range(config['interval_split'])))
    params:
        outpref=lambda wildcards: tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/1-23.processed.split'.format(**wildcards)
    run:
        import platform
        if platform.system() == 'Darwin':
            split_cmd = 'gsplit'
        elif platform.system() == 'Linux':
            split_cmd = 'split'
        else:
            assert(True, 'Error: platform must be Darwin or Linux')

        shell(split_cmd + ' -a 3 --additional-suffix=.tsv -d -n l/128 {input} {params.outpref}')

rule split_rezip:
    ''' Re zip the split file
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{version}/{cell}/1-23.processed.split{i}.tsv'
    output:
        config['out_dir'] + '/interval/pchic/javierre2016/{version}/{cell}/1-23.processed.split{i}.tsv.gz'
    shell:
        'gzip -c {input} > {output}'
