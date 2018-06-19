rule andersson2014_download:
    ''' Retrieves Andersson et al 2014 Fantom5
    '''
    input:
        HTTP.remote('http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed')
    output:
        tmpdir + '/interval/fantom5/andersson2014/unspecified/enhancer_tss_associations.bed'
    shell:
        'cp {input} {output}'

rule andersson2014_to_final:
    ''' Produces finalised output of fantom5 intervals
    '''
    input:
        bed = tmpdir + '/interval/fantom5/andersson2014/unspecified/enhancer_tss_associations.bed',
        gtf = tmpdir + '/Homo_sapiens.GRCh37.87.gtf.gz'
    output:
        tmpdir + '/interval/fantom5/andersson2014/unspecified/1-23.processed.tsv.gz'
    shell:
        'python scripts/andersson2014_to_final.py '
        '--inf {input.bed} '
        '--outf {output} '
        '--gtf {input.gtf} '
        '--cell_name Unspecified'

rule andersson2014_final_to_gcs:
    ''' Copy to google cloud storage
    '''
    input:
        tmpdir + '/interval/fantom5/andersson2014/unspecified/1-23.processed.tsv.gz'
    output:
        GS.remote(
            '{bucket}/{gs_dir}/interval/fantom5/andersson2014/unspecified/1-23.processed.tsv.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir']))
    shell:
        'cp {input} {output}'

rule andersson2014_raw_to_gcs:
    ''' Copy to google cloud storage
    '''
    input:
        tmpdir + '/interval/fantom5/andersson2014/unspecified/enhancer_tss_associations.bed'
    output:
        GS.remote(
            '{bucket}/{gs_dir}/interval/fantom5/andersson2014/unspecified/1-23.raw.bed.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir']))
    shell:
        'cp {input} {output}'
