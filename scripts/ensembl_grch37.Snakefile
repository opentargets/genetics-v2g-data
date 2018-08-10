from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

rule ensembl_download_gtf:
    ''' Downloads the Ensembl GRCh37 GTF
    '''
    input:
        FTPRemoteProvider().remote('ftp.ensembl.org/pub/grch37/update/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz',
                                   keep_local=False)
    output:
        tmpdir + '/Homo_sapiens.GRCh37.87.gtf.gz'
    shell:
        'cp {input} {output}'

rule ensembl_gtf_to_tssbed:
    ''' Converts gtf to TSS bed file
    Notes:
        - Keeps protein coding transcripts only.
    '''
    input:
        tmpdir + '/Homo_sapiens.GRCh37.87.gtf.gz'
    output:
        tmpdir + '/Homo_sapiens.GRCh37.87.tss.bed.gz'
    shell:
        'python scripts/ensembl_gtf_to_tss_bed.py '
        '--inf {input} '
        '--outf {output}'

rule ensembl_to_uniprot_download:
    ''' Uses biomart to download file to map from Ensembl ID to Uniprot ID
    '''
    output:
        tmpdir + '/Ensembl_to_hgnc_map.GRCh37.tsv.gz'
    shell:
        'wget -O - \'http://grch37.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "transcription_start_site" /><Attribute name = "hgnc_symbol" /><Attribute name = "uniprot_gn" /><Attribute name = "chromosome_name" /></Dataset></Query>\' | gzip -c > {output}'
