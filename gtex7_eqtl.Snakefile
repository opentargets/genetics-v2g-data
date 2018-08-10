#!/usr/bin/env snakemake
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
import pandas as pd
from pprint import pprint
from datetime import date

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']
version = date.today().strftime("%y%m%d")

## Make targets for GTEX V7 eQTL dataset
targets = []
tissues, = GSRemoteProvider().glob_wildcards(config['gtex7_raw_gs_dir'] + '/{tissues}.v7.signif_variant_gene_pairs.txt.gz')

# Create cis-regulatory data target
for tissue in tissues:
    target = '{out_dir}/{data_type}/{exp_type}/{source}/{version}/{cell_type}/{chrom}.{proc}.processed.tsv.gz'.format(
        out_dir=config['out_dir'],
        data_type='qtl',
        exp_type='eqtl',
        source='gtex_v7',
        version=version,
        cell_type=tissue,
        chrom='1-23',
        proc='cis_reg')
    targets.append(target)

# "all" must be first rule that is encounterd
rule all:
    ''' Master rule to trigger all targets (defined above) '''
    input:
        targets

rule download_from_gcs:
    ''' Copies from gcs to tmp folder
    '''
    output:
        temp(tmpdir + '/qtl/gtex_download/{tissue}.v7.signif_variant_gene_pairs.txt.gz')
    params:
        input=lambda wildcards: config['gtex7_raw_gs_dir'] + '/{tissue}.v7.signif_variant_gene_pairs.txt.gz'.format(**wildcards)
    shell:
        'gsutil cp {params.input} {output}'

rule format_gtex:
    ''' Formats GTEx v7 significant_pairs file in standard way
    '''
    input:
        tmpdir + '/qtl/gtex_download/{tissue}.v7.signif_variant_gene_pairs.txt.gz'
    output:
        '{out_dir}/{data_type}/{exp_type}/{source}/{version}/{{tissue}}/{chrom}.{proc}.processed.tsv.gz'.format(
            out_dir=config['out_dir'],
            data_type='qtl',
            exp_type='eqtl',
            source='gtex_v7',
            version=version,
            chrom='1-23',
            proc='cis_reg')
    shell:
        'python scripts/gtex7_format_cisreg.py '
        '--inf {input} '
        '--outf {output}'
