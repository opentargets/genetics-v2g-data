rule thurman2012_download:
    ''' Retrieves DHS correlation data from the Ensembl ENCODE FTP
    '''
    input:
        HTTP.remote('http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/dhs_gene_connectivity/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz')
    output:
        tmpdir + '/interval/dhscor/thurman2012/unspecified/dhs_correlations.bed.gz'
    shell:
        'cp {input} {output}'

rule thurman2012_to_final:
    ''' Produces finalised output of fantom5 intervals
    '''
    input:
        bed = tmpdir + '/interval/dhscor/thurman2012/unspecified/dhs_correlations.bed.gz',
        gtf = tmpdir + '/Homo_sapiens.GRCh37.87.gtf.gz'
    output:
        tmpdir + '/interval/dhscor/thurman2012/unspecified/1-23.processed.tsv.gz'
    shell:
        'python scripts/thurman2012_to_final.py '
        '--inf {input.bed} '
        '--outf {output} '
        '--gtf {input.gtf} '
        '--cell_name Unspecified'

rule thurman2012_final_to_gcs:
    ''' Copy final to google cloud storage
    '''
    input:
        tmpdir + '/interval/dhscor/thurman2012/unspecified/1-23.processed.tsv.gz'
    output:
        GS.remote(
            '{bucket}/{gs_dir}/interval/dhscor/thurman2012/unspecified/1-23.processed.tsv.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir']))
    shell:
        'cp {input} {output}'

rule thurman2012_raw_to_gcs:
    ''' Copy raw to google cloud storage
    '''
    input:
        tmpdir + '/interval/dhscor/thurman2012/unspecified/dhs_correlations.bed.gz'
    output:
        GS.remote(
            '{bucket}/{gs_dir}/interval/dhscor/thurman2012/unspecified/1-23.raw.bed.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir']))
    shell:
        'cp {input} {output}'
