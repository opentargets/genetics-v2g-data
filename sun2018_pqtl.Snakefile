#!/usr/bin/env snakemake
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import pandas as pd
from pprint import pprint
import gzip
from datetime import date

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']
version = date.today().strftime("%y%m%d")
KEEP_LOCAL = False

# Initiate remote handles
GS = GSRemoteProvider()

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

# DEBUG
# sun2018_manifest = sun2018_manifest.loc[:2, :]

# Get list of ids that exist on GCS
# gcs_ids, files = GS.glob_wildcards(config['sun2018_in_gs_dir'] + '/UBERON_0001969/{uniprot}/{file}.tsv.gz')
# gcs_ids = set(gcs_ids)
# to_keep = sun2018_manifest['gcs_id'].isin(gcs_ids)
# print(sun2018_manifest.shape)
# sun2018_manifest = sun2018_manifest.loc[to_keep, :]
# print(sun2018_manifest.shape)

# Create cis-regulatory data target
target = '{out_dir}/{data_type}/{exp_type}/{source}/{version}/{cell_type}/{chrom}.pval{pval}.{proc}.processed.tsv.gz'.format(
    out_dir=config['out_dir'],
    data_type='qtl',
    exp_type='pqtl',
    source='sun2018',
    version=version,
    cell_type='Blood_plasma',
    chrom='1-22',
    proc='cis_reg',
    pval=config['sun2018_cis_pval'])
targets.append(target)

# Create trans-regulatory data target
# target = '{bucket}/{gs_dir}/{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.pval{pval}.{proc}.processed.tsv.gz'.format(
#     bucket=config['gs_bucket'],
#     gs_dir=config['gs_dir'],
#     data_type='qtl',
#     exp_type='pqtl',
#     source='sun2018',
#     cell_type='Blood_plasma',
#     chrom='1-22',
#     proc='trans_reg',
#     pval=config['sun2018_trans_pval'])
# if UPLOAD:
#     targets.append(GS.remote(target))

# "all" must be first rule that is encounterd
rule all:
    ''' Master rule to trigger all targets (defined above) '''
    input:
        targets

# Import workflows
# include: 'scripts/ensembl_grch37.Snakefile'

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
        # lambda wildcards: GSRemoteProvider().remote('{gs_dir}/UBERON_0001969/{uniprot_id}/{chrom}-SUN2018-UBERON_0001969-{uniprot_id}.tsv.gz'.format(
        #         gs_dir=config['sun2018_in_gs_dir'],
        #         uniprot_id=get_manifest_info(sun2018_manifest, wildcards.ensembl_id,
        #                                     'gcs_id'),
        #         chrom=get_manifest_info(sun2018_manifest, wildcards.ensembl_id,
        #                                 'chrom')), keep_local=KEEP_LOCAL)
    output:
        temp('{tmpdir}/qtl/pqtl/sun2018/{{version}}/{{cell_type}}/{{ensembl_id}}.pval{pval}.cis_reg.tsv.gz'.format(
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

rule concate_dataset:
    ''' Concatenate separate gene files into a single file, and upload to GCS
    '''
    input:
        ['{tmpdir}/qtl/pqtl/sun2018/{{version}}/{{cell_type}}/{ensembl_id}.pval{{pval}}.{{proc}}.tsv.gz'.format(
                tmpdir=tmpdir,
                ensembl_id=ensembl_id)
         for ensembl_id in list(sun2018_manifest['ensembl_id'].unique())]
    output:
        '{out_dir}/{{data_type}}/{{exp_type}}/{{source}}/{{version}}/{{cell_type}}/{{chrom}}.pval{{pval}}.{{proc}}.processed.tsv.gz'.format(
            out_dir=config['out_dir'])
    run:
        with gzip.open(output[0], 'w') as out_h:
            for i, inf in enumerate(input):
                # Get ensembl ID name from path
                ensembl_id = os.path.split(inf)[1].split('.')[0]
                # Open file
                with gzip.open(inf, 'r') as in_h:
                    # Write header if first file
                    header = in_h.readline().decode().rstrip().split('\t')
                    if i == 0:
                        header.insert(4, 'ensembl_id')
                        out_h.write(('\t'.join(header) + '\n').encode())
                    # Write all lines to file
                    for line in in_h:
                        line = line.decode().rstrip().split('\t')
                        line.insert(4, ensembl_id)
                        out_h.write(('\t'.join(line) + '\n').encode())

# rule extract_trans_data:
#     ''' Extracts trans-regulatory data from Sun pQTL data
#     '''
#     input:
#         lambda wildcards: [GSRemoteProvider().remote('{bucket}/{gs_dir}/{soma_id}/{soma_id}_chrom_{chrom}_meta_final_v1.tsv.gz'.format(
#                 bucket=config['gs_bucket'],
#                 gs_dir=config['sun2018_raw_gs_dir'],
#                 soma_id=get_manifest_info(sun2018_manifest, wildcards.ensembl_id, 'SOMAMER_ID'),
#                 chrom=str(chrom))) for chrom in range(1, 23)]
#     output:
#         temp('{tmpdir}/qtl/pqtl/sun2018/{{cell_type}}/{{ensembl_id}}.pval{pval}.trans_reg.tsv.gz'.format(
#                 tmpdir=tmpdir,
#                 pval=config['sun2018_trans_pval']))
#     params:
#         chrom=lambda wildcards: get_manifest_info(sun2018_manifest, wildcards.ensembl_id, 'chrom'),
#         tss=lambda wildcards: int(float(get_manifest_info(sun2018_manifest, wildcards.ensembl_id, 'tss'))),
#         window=int(config['sun2018_cis_window']),
#         pval=float(config['sun2018_trans_pval'])
#     shell:
#         'python scripts/sun2018_extract_trans_data.py '
#         '--infs {input} '
#         '--outf {output} '
#         '--chrom {params.chrom} '
#         '--tss_pos {params.tss} '
#         '--cis_window {params.window} '
#         '--pval {params.pval}'
