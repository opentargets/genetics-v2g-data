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

## Make targets for GTEX V7 eQTL dataset

Load tisues from manifest
tissues, = GS.glob_wildcards(config['gtex7_raw_gs_dir'] + '/{tissues}.v7.signif_variant_gene_pairs.txt.gz')

# Create cis-regulatory data target
for tissue in tissues:
    target = '{bucket}/{gs_dir}/{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.{proc}.processed.tsv.gz'.format(
        bucket=config['gs_bucket'],
        gs_dir=config['gs_dir'],
        data_type='qtl',
        exp_type='eqtl',
        source='gtex_v7',
        cell_type=tissue,
        chrom='1-23',
        proc='cis_reg')
   if UPLOAD:
       targets.append(GS.remote(target))

# "all" must be first rule that is encounterd
rule all:
    ''' Master rule to trigger all targets (defined above) '''
    input:
        targets

# Import workflows
include: 'scripts/ensembl_grch37.Snakefile'

rule format_gtex:
    ''' Formats GTEx v7 significant_pairs file in standard way
    '''
    input:
        GSRemoteProvider().remote(
            '{gs_dir}/{{tissue}}.v7.signif_variant_gene_pairs.txt.gz'.format(
                gs_dir=config['gtex7_raw_gs_dir']))
    output:
        GSRemoteProvider().remote(
            '{bucket}/{gs_dir}/qtl/eqtl/gtex_v7/{{tissue}}/1-23.cis_reg.processed.tsv.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir']))
    shell:
        'python scripts/gtex7_format_cisreg.py '
        '--inf {input} '
        '--outf {output}'
