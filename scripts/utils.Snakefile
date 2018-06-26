rule unzip_processed:
    ''' Unzips a processed file and removes the header
    '''
    input:
        GSRemoteProvider().remote('{filename}.processed.tsv.gz')
    output:
        temp('{filename}.processed.tsv')
    shell:
        'gunzip -c {input} | tail -n +2 > {output}'

rule split_processed:
    ''' Splits an unzipped file
    '''
    input:
        '{filename}.processed.tsv'
    output:
        [ GSRemoteProvider().remote('{{filename}}.processed.split{i:03d}.tsv.gz'.format(i=i))
          for i in range(config['interval_split']) ]
    params:
        outpref= lambda wildcards: '{f}.processed.split'.format(f=wildcards.filename),
        chunks=config['interval_split']
    shell:
        "gsplit -a 3 --additional-suffix=.tsv --filter='gzip > $FILE.gz' -d -n l/{params.chunks} {input} {params.outpref}"
