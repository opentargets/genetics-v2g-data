#!/usr/bin/env snakemake
'''
Ed Mountjoy (June 2018)

Master script to create targets for all datasets. Workflow, based of POSTGAP
makefile: https://github.com/Ensembl/postgap/blob/master/Makefile

'''

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import pandas as pd
from pprint import pprint

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']

# Initiate remote handles
GS = GSRemoteProvider()
# FTP = FTPRemoteProvider()
# HTTP = HTTPRemoteProvider()

targets = []

#
# Make QTL dataset targets -----------------------------------------------------
#

### Make targets for GTEX V7 eQTL dataset

# Load tisues from manifest
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
    # targets.append(GS.remote(target))

### Make targets for Sun et al pQTL dataset

# Load manifest
sun2018_manifest = pd.read_csv(config['sun2018_manifest'], sep='\t', header=0)
# Drop non automsomal genes
valid_chroms = set([str(chrom) for chrom in range(1, 23)])
sun2018_manifest.chrom = sun2018_manifest.chrom.astype(str)
sun2018_manifest = sun2018_manifest.loc[sun2018_manifest.chrom.isin(valid_chroms), :]
# Drop duplicate genes
sun2018_manifest = sun2018_manifest.drop_duplicates(subset='gene')


# Create cis-regulatory data target
target = '{bucket}/{gs_dir}/{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.pval{pval}.{proc}.processed.tsv.gz'.format(
    bucket=config['gs_bucket'],
    gs_dir=config['gs_dir'],
    data_type='qtl',
    exp_type='pqtl',
    source='sun2018',
    cell_type='Blood_plasma',
    chrom='1-22',
    proc='cis_reg',
    pval=config['sun2018_cis_pval'])
# targets.append(GS.remote(target))
targets.append(target)

# Create trans-regulatory data target
target = '{bucket}/{gs_dir}/{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.pval{pval}.{proc}.processed.tsv.gz'.format(
    bucket=config['gs_bucket'],
    gs_dir=config['gs_dir'],
    data_type='qtl',
    exp_type='pqtl',
    source='sun2018',
    cell_type='Blood_plasma',
    chrom='1-22',
    proc='trans_reg',
    pval=config['sun2018_trans_pval'])
# targets.append(GS.remote(target))

#
# Make interval dataset targets ------------------------------------------------
#

## Make targets for Javierre 2016 PCHiC

# Get list of cell line names
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
    # targets.append(GS.remote(target))

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
        # targets.append(GS.remote(target))

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
    # targets.append(GS.remote(target))

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
# targets.append(GS.remote(target))

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
    # targets.append(GS.remote(target))

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
# targets.append(GS.remote(target))

### Make targets for Thurman 2012 DHS-promoter correlation dataset

# Processed
target = '{bucket}/{gs_dir}/{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.{proc}.tsv.gz'.format(
    bucket=config['gs_bucket'],
    gs_dir=config['gs_dir'],
    data_type='interval',
    exp_type='dhscor',
    source='thurman2012',
    cell_type='unspecified',
    proc='processed',
    chrom='1-23')
# targets.append(GS.remote(target))

# Split files
for i in range(config['interval_split']):
    target = '{bucket}/{gs_dir}/{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.{proc}.split{i:03d}.tsv.gz'.format(
        bucket=config['gs_bucket'],
        gs_dir=config['gs_dir'],
        data_type='interval',
        exp_type='dhscor',
        source='thurman2012',
        cell_type='unspecified',
        proc='processed',
        chrom='1-23',
        i=i)
    # targets.append(GS.remote(target))

# Raw
target = '{bucket}/{gs_dir}/{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.{proc}.bed.gz'.format(
    bucket=config['gs_bucket'],
    gs_dir=config['gs_dir'],
    data_type='interval',
    exp_type='dhscor',
    source='thurman2012',
    cell_type='unspecified',
    proc='raw',
    chrom='1-23')
# targets.append(GS.remote(target))

pprint(targets)

#
# Run rules to produce targets -------------------------------------------------
#

# "all" must be first rule that is encounterd
rule all:
    ''' Master rule to trigger all targets (defined above) '''
    input:
        targets

# Import workflows
include: 'scripts/andersson2014_fantom5.Snakefile'
include: 'scripts/ensembl_grch37.Snakefile'
include: 'scripts/javierre2016_pchic.Snakefile'
include: 'scripts/sun2018_pqtl.Snakefile'
include: 'scripts/thurman2012_dhscor.Snakefile'
include: 'scripts/gtex7_eqtl.Snakefile'
include: 'scripts/utils.Snakefile'
