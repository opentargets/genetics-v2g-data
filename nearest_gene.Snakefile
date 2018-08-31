#!/usr/bin/env snakemake
import pandas as pd
from pprint import pprint
from datetime import date
import gzip

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']
version = date.today().strftime("%y%m%d")

# Make targets
targets = []

# Processed
target = '{out_dir}/nearest_gene/{version}/nearest_gene.tsv.gz'.format(
    out_dir=config['out_dir'],
    version=version)
targets.append(target)

# "all" must be first rule that is encounterd
rule all:
    ''' Master rule to trigger all targets (defined above) '''
    input:
        targets

# Import workflows
include: 'scripts/ensembl_grch37.Snakefile'

rule download_from_gcs:
    ''' Copies from gcs to tmp folder
    '''
    output:
        temp(tmpdir + '/nearest_gene/{version}/homo_sapiens_incl_consequences.vcf.gz')
    params:
        input='gs://genetics-portal-data/homo_sapiens_incl_consequences.vcf.gz'
    shell:
        'gsutil cp {params.input} {output}'

rule make_variant_index_bed:
    ''' Converts VCF to bed format. Warning this is slow ~2 hour!
    '''
    input:
        tmpdir + '/nearest_gene/{version}/homo_sapiens_incl_consequences.vcf.gz'
    output:
        tmpdir + '/nearest_gene/{version}/homo_sapiens_incl_consequences.bed.gz'
    shell:
        'pypy3 scripts/vcf_to_bed.py --invcf {input} --outbed {output}'

rule find_nearest_protein_coding:
    ''' Finds nearest protein coding gene for each variant
    '''
    input:
        vars=tmpdir + '/nearest_gene/{version}/homo_sapiens_incl_consequences.bed.gz',
        genes=tmpdir + '/Homo_sapiens.GRCh37.87.tss.protein_coding.bed.gz'
    output:
        tmpdir + '/nearest_gene/{version}/nearest_gene.protein_coding.tsv.gz'
    resources:
        threads=2
    shell:
        'bedtools closest -d -t first -a {input.vars} -b {input.genes} '
        '| cut -f 4,9,11 | gzip -c > {output}'

rule find_nearest_any:
    ''' Finds nearest (any) gene for each variant
    '''
    input:
        vars=tmpdir + '/nearest_gene/{version}/homo_sapiens_incl_consequences.bed.gz',
        genes=tmpdir + '/Homo_sapiens.GRCh37.87.tss.bed.gz'
    output:
        tmpdir + '/nearest_gene/{version}/nearest_gene.any.tsv.gz'
    resources:
        threads=2
    shell:
        'bedtools closest -d -t first -a {input.vars} -b {input.genes} '
        '| cut -f 4,9,11 | gzip -c > {output}'

rule merge_nearest_proteincoding_and_any:
    ''' Merges the protein coding results with the any gene results
    '''
    input:
        protein_coding=tmpdir + '/nearest_gene/{version}/nearest_gene.protein_coding.tsv.gz',
        any=tmpdir + '/nearest_gene/{version}/nearest_gene.any.tsv.gz'
    output:
        '{out_dir}/nearest_gene/{{version}}/nearest_gene.tsv.gz'.format(out_dir=config['out_dir'])
    run:
        # Open all handles
        with gzip.open(input['protein_coding'], 'rt') as in_pc_h, \
             gzip.open(input['any'], 'rt') as in_any_h, \
             gzip.open(output[0], 'wt') as out_h:

            # Write header
            header = ['varid',
                      'ensemblid_protein_coding',
                      'distance_protein_coding',
                      'ensemblid_any',
                      'distance_any']
            out_h.write('\t'.join(header) + '\n')

            # Iterate over input file
            for pc_line in in_pc_h:
                # Parse data
                pc_var, pc_gene, pc_dist = pc_line.rstrip().split('\t')
                any_var, any_gene, any_dist = in_any_h.readline().rstrip().split('\t')
                # Check that varids are the same
                if not pc_var == any_var:
                    sys.exit('Error: files are not in sync (pc_var != any_var)')
                # Write lines
                out_row = [pc_var, pc_gene, pc_dist, any_gene, any_dist]
                out_h.write('\t'.join(out_row) + '\n')
