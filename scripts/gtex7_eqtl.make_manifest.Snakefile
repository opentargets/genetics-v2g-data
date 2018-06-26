# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']

rule all:
    input:
        config['gtex7_manifest']

rule make_manifest:
    ''' Downloads tissues names to make a manifest for GTEx dataset
    '''
    output:
        config['gtex7_manifest']
    params:
        gs_dir=config['gtex7_raw_gs_dir']
    shell:
        "gsutil ls 'gs://{params.gs_dir}/*.v7.signif_variant_gene_pairs.txt.gz' | cut -f 5 -d '/' | sed 's/.v7.signif_variant_gene_pairs.txt.gz//g' > {output}"
