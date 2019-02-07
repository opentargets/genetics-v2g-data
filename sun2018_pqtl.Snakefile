#!/usr/bin/env snakemake
import pandas as pd
from pprint import pprint
import gzip
from datetime import date

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']
version = date.today().strftime("%y%m%d")
KEEP_LOCAL = False

targets = []

## Make targets for Sun et al pQTL dataset

# Load manifest
sun2018_manifest = pd.read_csv(config['sun2018_manifest_local'], sep='\t', header=0)
# Drop non automsomal genes
valid_chroms = set([str(chrom) for chrom in range(1, 23)])
sun2018_manifest.chrom = sun2018_manifest.chrom.astype(str)
sun2018_manifest = sun2018_manifest.loc[sun2018_manifest.chrom.isin(valid_chroms), :]
# Drop uniprots that don't map to an Ensmebl
sun2018_manifest = sun2018_manifest.dropna(subset=['gene'])
# Keep first if Ensembl gene maps to multiple rows
sun2018_manifest = sun2018_manifest.drop_duplicates(subset='gene')
# Make identifier as it appears on GCS
sun2018_manifest['gcs_id'] = sun2018_manifest.UniProt.apply(lambda x: 'UNIPROT_' + '_'.join(x.split(',')))

# Create cis-regulatory data target
target = '{out_dir}/{data_type}/{exp_type}/{source}/{version}/data.parquet'.format(
    out_dir=config['out_dir'],
    data_type='qtl',
    exp_type='pqtl',
    source='sun2018',
    version=version,
    cell_type=config['sun2018_cell'])
targets.append(target)

# "all" must be first rule that is encounterd
rule all:
    ''' Master rule to trigger all targets (defined above) '''
    input:
        targets

def get_manifest_info(manifest, ensembl_id, key):
    ''' Extract info from manifest for a given ensembl_id
    Args:
        manifest (DF): Sun et al manifest file
        ensembl_id (str): ensembl id
        key (str): column to return
    Returns:
        value (str)
    '''
    row = manifest.loc[manifest.ensembl_id == ensembl_id, :].iloc[0, :]
    return row[key]

rule download_from_gcs:
    ''' Copies from gcs to tmp folder
    '''
    output:
        temp(tmpdir + '/qtl/pqtl_download/UBERON_0001969/{uniprot_id}/{chrom}-SUN2018-UBERON_0001969-{uniprot_id}.tsv.gz')
    params:
        input=lambda wildcards: config['sun2018_in_gs_dir'] + '/UBERON_0001969/{uniprot_id}/{chrom}-SUN2018-UBERON_0001969-{uniprot_id}.tsv.gz'.format(**wildcards)
    shell:
        'gsutil cp {params.input} {output}'

rule extract_cis_data:
    ''' Extracts cis-regulatory data from Sun pQTL data for a single gene
    '''
    input:
        lambda wildcards: tmpdir + '/qtl/pqtl_download/UBERON_0001969/{uniprot_id}/{chrom}-SUN2018-UBERON_0001969-{uniprot_id}.tsv.gz'.format(
                gs_dir=config['sun2018_in_gs_dir'],
                uniprot_id=get_manifest_info(sun2018_manifest, wildcards.ensembl_id, 'gcs_id'),
                chrom=get_manifest_info(sun2018_manifest, wildcards.ensembl_id, 'chrom'))
    output:
        temp('{tmpdir}/qtl/pqtl/sun2018/{{version}}/{{ensembl_id}}.cis.tsv.gz'.format(
                tmpdir=tmpdir,
                pval=config['sun2018_cis_pval']))
    params:
        chrom=lambda wildcards: get_manifest_info(sun2018_manifest, wildcards.ensembl_id, 'chrom'),
        tss=lambda wildcards: int(float(get_manifest_info(sun2018_manifest, wildcards.ensembl_id, 'tss'))),
        window=int(config['sun2018_cis_window']),
        pval=float(config['sun2018_cis_pval'])
    shell:
        'python scripts/sun2018_extract_cis_data.py '
        '--inf {input} '
        '--outf {output} '
        '--chrom {params.chrom} '
        '--tss_pos {params.tss} '
        '--cis_window {params.window} '
        '--pval {params.pval}'

rule to_parquet:
    ''' Concatenate separate gene files into a single file, and upload to GCS
    '''
    input:
        data = ['{tmpdir}/qtl/pqtl/sun2018/{{version}}/{ensembl_id}.cis.tsv.gz'.format(
                tmpdir=tmpdir,
                ensembl_id=ensembl_id)
         for ensembl_id in list(sun2018_manifest['ensembl_id'].unique()) ]
         # for ensembl_id in ['ENSG00000160712']], # DEBUG
    output:
        directory('{out_dir}/{{data_type}}/{{exp_type}}/{{source}}/{{version}}/data.parquet'.format(
            out_dir=config['out_dir']))
    params:
        cell = config['sun2018_cell'],
        pattern = '{tmpdir}/qtl/pqtl/sun2018/{version}/*.cis.tsv.gz'.format(
            tmpdir=tmpdir,
            version=version)
    shell:
        'python scripts/sun2018_to_parquet.py '
        '--infs {params.pattern} '
        '--outf {output} '
        '--cell_code {params.cell}'
