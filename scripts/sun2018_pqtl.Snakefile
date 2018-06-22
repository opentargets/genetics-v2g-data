import gzip

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

rule extract_cis_data:
    ''' Extracts cis-regulatory data from Sun pQTL data for a single gene
    '''
    input:
        lambda wildcards: GSRemoteProvider().remote(
            '{bucket}/{gs_dir}/{soma_id}/{soma_id}_chrom_{chrom}_meta_final_v1.tsv.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['sun2018_raw_gs_dir'],
                soma_id=get_manifest_info(sun2018_manifest, wildcards.ensembl_id,
                                          'SOMAMER_ID'),
                chrom=get_manifest_info(sun2018_manifest, wildcards.ensembl_id,
                                        'chrom')))
    output:
        temp('{tmpdir}/qtl/pqtl/sun2018/{{cell_type}}/{{ensembl_id}}.pval{pval}.cis_reg.tsv.gz'.format(
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

rule extract_trans_data:
    ''' Extracts trans-regulatory data from Sun pQTL data
    '''
    input:
        lambda wildcards: [GSRemoteProvider().remote('{bucket}/{gs_dir}/{soma_id}/{soma_id}_chrom_{chrom}_meta_final_v1.tsv.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['sun2018_raw_gs_dir'],
                soma_id=get_manifest_info(sun2018_manifest, wildcards.ensembl_id, 'SOMAMER_ID'),
                chrom=str(chrom))) for chrom in range(1, 23)]
    output:
        temp('{tmpdir}/qtl/pqtl/sun2018/{{cell_type}}/{{ensembl_id}}.pval{pval}.trans_reg.tsv.gz'.format(
                tmpdir=tmpdir,
                pval=config['sun2018_trans_pval']))
    params:
        chrom=lambda wildcards: get_manifest_info(sun2018_manifest, wildcards.ensembl_id, 'chrom'),
        tss=lambda wildcards: int(float(get_manifest_info(sun2018_manifest, wildcards.ensembl_id, 'tss'))),
        window=int(config['sun2018_cis_window']),
        pval=float(config['sun2018_trans_pval'])
    shell:
        'python scripts/sun2018_extract_trans_data.py '
        '--infs {input} '
        '--outf {output} '
        '--chrom {params.chrom} '
        '--tss_pos {params.tss} '
        '--cis_window {params.window} '
        '--pval {params.pval}'

rule concate_dataset:
    ''' Concatenate separate gene files into a single file, and upload to GCS
    '''
    input:
        ['{tmpdir}/qtl/pqtl/sun2018/{{cell_type}}/{ensembl_id}.pval{{pval}}.{{proc}}.tsv.gz'.format(
                tmpdir=tmpdir,
                # pval=config['sun2018_cis_pval'],
                ensembl_id=ensembl_id)
         for ensembl_id in list(sun2018_manifest['ensembl_id'].unique())]
    output:
        GSRemoteProvider().remote('{bucket}/{gs_dir}/qtl/pqtl/sun2018/{{cell_type}}/{{chrom}}.pval{{pval}}.{{proc}}.processed.tsv.gz'.format(
            bucket=config['gs_bucket'],
            gs_dir=config['gs_dir']))
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
