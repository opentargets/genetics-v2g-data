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
