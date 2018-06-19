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
    ''' Extracts cis-regulatory data from Sun pQTL data
    '''
    input:
        lambda wildcards: GS.remote(
            '{bucket}/{gs_dir}/{soma_id}/{soma_id}_chrom_{chrom}_meta_final_v1.tsv.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['sun2018_raw_gs_dir'],
                soma_id=get_manifest_info(sun2018_manifest,
                                          wildcards.ensembl_id,
                                          'SOMAMER_ID'),
                chrom=wildcards.chrom))
    output:
        GS.remote('{bucket}/{gs_dir}/qtl/pqtl/sun2018/{{cell_type}}/{{ensembl_id}}/{{chrom}}.pval{pval}.cis_reg.tsv.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir'],
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
        lambda wildcards: [GS.remote('{bucket}/{gs_dir}/{soma_id}/{soma_id}_chrom_{chrom}_meta_final_v1.tsv.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['sun2018_raw_gs_dir'],
                soma_id=get_manifest_info(sun2018_manifest, wildcards.ensembl_id, 'SOMAMER_ID'),
                chrom=str(chrom))) for chrom in range(1, 23)]
    output:
        GS.remote('{bucket}/{gs_dir}/qtl/pqtl/sun2018/{{cell_type}}/{{ensembl_id}}/{{chrom}}.pval{pval}.trans_reg.tsv.gz'.format(
            bucket=config['gs_bucket'],
            gs_dir=config['gs_dir'],
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

rule format_full_data:
    ''' Reformats the raw data file to a standard format
    '''
    input:
        lambda wildcards: GS.remote(
            '{bucket}/{gs_dir}/{soma_id}/{soma_id}_chrom_{chrom}_meta_final_v1.tsv.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['sun2018_raw_gs_dir'],
                soma_id=get_manifest_info(sun2018_manifest,
                                          wildcards.ensembl_id,
                                          'SOMAMER_ID'),
                chrom=wildcards.chrom))
    output:
        GS.remote('{bucket}/{gs_dir}/qtl/pqtl/sun2018/{{cell_type}}/{{ensembl_id}}/{{chrom}}.pval1.all.tsv.gz'.format(
            bucket=config['gs_bucket'],
            gs_dir=config['gs_dir']))
    shell:
        'python scripts/sun2018_format_full_data.py '
        '--inf {input} '
        '--outf {output}'
