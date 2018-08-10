#!/usr/bin/env snakemake
# from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
# from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
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
target = '{out_dir}/{data_type}/{exp_type}/{source}/{version}/{cell_type}/{chrom}.{proc}.tsv.gz'.format(
    out_dir=config['out_dir'],
    data_type='interval',
    exp_type='fantom5',
    source='andersson2014',
    version=version,
    cell_type='unspecified',
    proc='processed',
    chrom='1-23')
targets.append(target)

# Split files
for i in range(config['interval_split']):
    target = '{out_dir}/{data_type}/{exp_type}/{source}/{version}/{cell_type}/{chrom}.{proc}.split{i:03d}.tsv.gz'.format(
        out_dir=config['out_dir'],
        data_type='interval',
        exp_type='fantom5',
        source='andersson2014',
        version=version,
        cell_type='unspecified',
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

rule andersson2014_download:
    ''' Retrieves Andersson et al 2014 Fantom5
    '''
    input:
        HTTPRemoteProvider().remote('http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed',
                                    keep_local=False)
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
        config['out_dir'] + '/interval/fantom5/andersson2014/{version}/unspecified/1-23.processed.tsv.gz'
    shell:
        'python scripts/andersson2014_to_final.py '
        '--inf {input.bed} '
        '--outf {output} '
        '--gtf {input.gtf} '
        '--cell_name Unspecified'

rule split_unzip_final:
    ''' Unzip and remove header in preparation for splitting
    '''
    input:
        config['out_dir'] + '/interval/fantom5/andersson2014/{version}/unspecified/1-23.processed.tsv.gz'
    output:
        temp(tmpdir + '/interval/fantom5/andersson2014/{version}/unspecified/1-23.processed.tsv')
    shell:
        'zcat < {input} | tail -n +2 > {output}'

rule split:
    ''' Splits file into many parts
    '''
    input:
        tmpdir + '/interval/fantom5/andersson2014/{version}/unspecified/1-23.processed.tsv'
    output:
        temp(expand(tmpdir + '/interval/fantom5/andersson2014/{{version}}/unspecified/1-23.processed.split{i:03d}.tsv',
               i=range(config['interval_split'])))
    params:
        outpref=lambda wildcards: tmpdir + '/interval/fantom5/andersson2014/{version}/unspecified/1-23.processed.split'.format(**wildcards)
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
        tmpdir + '/interval/fantom5/andersson2014/{version}/unspecified/1-23.processed.split{i}.tsv'
    output:
        config['out_dir'] + '/interval/fantom5/andersson2014/{version}/unspecified/1-23.processed.split{i}.tsv.gz'
    shell:
        'gzip -c {input} > {output}'
