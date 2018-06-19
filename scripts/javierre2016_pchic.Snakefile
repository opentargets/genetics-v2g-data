rule javierre2016_download:
    ''' Retrieves Javierre 2016 PCHiC file from the ebi FTP
    '''
    input:
        FTP.remote('ftp.ebi.ac.uk/pub/contrib/pchic/CHiCAGO/{cell}.merged_samples_12Apr2015_full.txt.gz')
    output:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.raw.gz'
    shell:
        'cp {input} {output}'

rule javierre2016_to_bed:
    ''' Formats Javierre file to bed.
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.raw.gz'
    output:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.raw.bed.gz'
    shell:
        'python scripts/javierre2016_to_bed.py '
        '--inf {input} '
        '--outf {output}'

rule javierre2016_tss_intersect:
    ''' Finds the intersect between PCHiC capture regions and gene TSS
    '''
    input:
        pchic = tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.raw.bed.gz',
        tss = tmpdir + '/Homo_sapiens.GRCh37.87.tss.bed.gz'
    output:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.tss.bed.gz'
    shell:
        'bedtools intersect -wa -wb '
        '-a {input.pchic} '
        '-b {input.tss} | '
        'gzip -c > {output}'

rule javierre2016_to_final:
    ''' Outputs finalised format for Javierre dataset
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.tss.bed.gz'
    output:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.final.tsv.gz'
    shell:
        'python scripts/javierre2016_to_final.py '
        '--inf {input} '
        '--outf {output} '
        '--cell_name {wildcards.cell}'

rule javierre2016_final_to_gcs:
    ''' Copy final to google cloud storage
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.final.tsv.gz'
    output:
        GS.remote(
            '{bucket}/{gs_dir}/interval/pchic/javierre2016/{{cell}}/1-23.processed.tsv.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir']))
    shell:
        'cp {input} {output}'

rule javierre2016_raw_to_gcs:
    ''' Copy raw to google cloud storage
    '''
    input:
        tmpdir + '/interval/pchic/javierre2016/{cell}/1-23.raw.gz'
    output:
        GS.remote(
            '{bucket}/{gs_dir}/interval/pchic/javierre2016/{{cell}}/1-23.raw.bed.gz'.format(
                bucket=config['gs_bucket'],
                gs_dir=config['gs_dir']))
    shell:
        'cp {input} {output}'
