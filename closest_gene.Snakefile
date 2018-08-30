#!/usr/bin/env snakemake
import pandas as pd
from pprint import pprint
from datetime import date

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']
version = date.today().strftime("%y%m%d")

# Make targets
targets = []

# Processed
target = '{out_dir}/closest_gene/{version}/closest_gene.tsv.gz'.format(
    out_dir=config['out_dir'],
    version=version)
targets.append(target)

# "all" must be first rule that is encounterd
rule all:
    ''' Master rule to trigger all targets (defined above) '''
    input:
        targets
        # 'tmp/Homo_sapiens.GRCh37.87.tss.protein_coding.bed.gz'
        # [tmpdir + '/closest_gene/{version}/homo_sapiens_incl_consequences.bed.gz'.format(version=version)]

# Import workflows
include: 'scripts/ensembl_grch37.Snakefile'

rule download_from_gcs:
    ''' Copies from gcs to tmp folder
    '''
    output:
        temp(tmpdir + '/closest_gene/{version}/homo_sapiens_incl_consequences.vcf.gz')
    params:
        input='gs://genetics-portal-data/homo_sapiens_incl_consequences.vcf.gz'
    shell:
        'gsutil cp {params.input} {output}'

rule make_variant_index_bed:
    ''' Converts VCF to bed format. Warning this is slow ~2 hour!
    '''
    input:
        tmpdir + '/closest_gene/{version}/homo_sapiens_incl_consequences.vcf.gz'
    output:
        tmpdir + '/closest_gene/{version}/homo_sapiens_incl_consequences.bed.gz'
    shell:
        'pypy3 scripts/vcf_to_bed.py --invcf {input} --outbed {output}'

rule find_closest_protein_coding:
    ''' Finds closest protein coding gene for each variant
    '''
    input:
        vars=tmpdir + '/closest_gene/{version}/homo_sapiens_incl_consequences.bed.gz',
        genes=tmpdir + '/Homo_sapiens.GRCh37.87.tss.protein_coding.bed.gz'
    output:
        tmpdir + '/closest_gene/{version}/closest_gene.protein_coding.tsv.gz'
    resources:
        threads=2
    shell:
        'bedtools closest -d -t first -a {input.vars} -b {input.genes} '
        '| cut -f 4,9,11 | gzip -c > {output}'

rule find_closest_any:
    ''' Finds closest (any) gene for each variant
    '''
    input:
        vars=tmpdir + '/closest_gene/{version}/homo_sapiens_incl_consequences.bed.gz',
        genes=tmpdir + '/Homo_sapiens.GRCh37.87.tss.bed.gz'
    output:
        tmpdir + '/closest_gene/{version}/closest_gene.any.tsv.gz'
    resources:
        threads=2
    shell:
        'bedtools closest -d -t first -a {input.vars} -b {input.genes} '
        '| cut -f 4,9,11 | gzip -c > {output}'

rule merge_closest_proteincoding_and_any:
    ''' Merges the protein coding results with the any gene results
    '''
    input:
        protein_coding=tmpdir + '/closest_gene/{version}/closest_gene.protein_coding.tsv.gz',
        any=tmpdir + '/closest_gene/{version}/closest_gene.any.tsv.gz'
    output:
        '{out_dir}/closest_gene/{{version}}/closest_gene.tsv.gz'.format(out_dir=config['out_dir'])
    run:
        # Load both dfs
        protein_coding = pd.read_csv(input['protein_coding'], sep='\t', header=None)
        protein_coding.columns = ['varid', 'ensemblid_protein_coding', 'distance_protein_coding']
        any_gene = pd.read_csv(input['any'], sep='\t', header=None)
        any_gene.columns = ['varid', 'ensemblid_any', 'distance_any']
        # Merge
        merged = pd.merge(protein_coding, any_gene, on='varid', how='outer')
        # Save
        merged.to_csv(output[0], sep='\t', index=None, compression='gzip')
