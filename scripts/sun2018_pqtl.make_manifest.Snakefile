from snakemake.remote.GS import RemoteProvider as GSRemoteProvider

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']

rule all:
    input:
        config['sun2018_manifest_local']

rule pqtl2018_download_manifest:
    ''' Sun et at pQTL manifest file from GS
    '''
    input:
        GSRemoteProvider().remote(config['sun2018_manifest_gs'], keep_local=False)
    output:
        temp(tmpdir + '/qtl/pqtl/sun2018/Sun_pQTL_SOMALOGIC_GWAS_protein_info.manifest.csv')
    shell:
        'cp {input} {output}'

rule pqtl2018_map_manifest:
    ''' Map Ensembl IDs to the uniprot IDs
        - In map file, caluclate average TSS position per gene
        - In manifest, split rows with multiple genes
        - Map uniprot to Ensembl
    '''
    input:
        mani=tmpdir + '/qtl/pqtl/sun2018/Sun_pQTL_SOMALOGIC_GWAS_protein_info.manifest.csv',
        map=tmpdir + '/Ensembl_to_hgnc_map.GRCh37.tsv.gz'
    output:
        config['sun2018_manifest_local']
    shell:
        'python scripts/sun2018_map_manifest.py --inf {input.mani} '
        '--map {input.map} '
        '--outf {output}'

rule ensembl_to_uniprot_download:
    ''' Uses biomart to download file to map from Ensembl ID to Uniprot ID
    '''
    output:
        tmpdir + '/Ensembl_to_hgnc_map.GRCh37.tsv.gz'
    shell:
        'wget -O - \'http://grch37.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "transcription_start_site" /><Attribute name = "hgnc_symbol" /><Attribute name = "uniprot_gn" /><Attribute name = "chromosome_name" /></Dataset></Query>\' | gzip -c > {output}'
