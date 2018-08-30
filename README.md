Get cis-regulatory datasets
===========================

Process cis-regulatory datasets and upload them to 'genetics-portal-data
' bucket.

#### Datasets

### Interval datasets
- PCHiC from Javierre et al 2016
  - [Publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5123897/)
  - [Dataset](ftp://ftp.ebi.ac.uk/pub/contrib/pchic/CHiCAGO/))
  - Score = CHiCAGO score (Cairns et al., 2016) (CHiCAGO >= 5 only)
- Fantom5 enhancer-TSS associations from Andersson et al 2014
  - [Publication](https://www.nature.com/articles/nature12787#linking-enhancer-usage-with-tss-expression)
  - [Dataset](http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed)
  - Score = R-squared (FDR > 0.05 only)
- DHS correlation with gene promoters from Thurman et al 2012
  - [Publication](https://www.nature.com/articles/nature11232#a-map-of-distal-dhstopromoter-connections)
  - [Dataset](http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/dhs_gene_connectivity/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz)
  - Score = R-squared (R-squared > 0.7 only)

Output file naming convention:
  `{data type}/{experiment type}/{source}/{cell|tissue type}/{chrom}.{processed|raw}.tsv.gz`

Interval output columns:
  - `chrom`, `start`, `end`: 1-based genomic co-ordinates. Start and end are inclusive.
  - `ensembl_id`: Ensembl gene ID
  - `score`: feature score from original dataset (unmodified)
  - `cell_type`: cell or tissue type of this experiment

### QTL datasets
- GTEx v7
  - Notes:
    - None
- pQTL from Sun et al 2018
  - TODO
  - [Publication](https://www.nature.com/articles/s41586-018-0175-2)
  - [Data available from](http://www.phpc.cam.ac.uk/ceu/proteins/)
  - Data downloaded to `gs://genetics-portal-data/pqtl_sun2018`
  - Notes:
    - Do we only keep cis-associations?
    - Do we need to exclude "pQTLs without evidence against binding effects"
      - Main body of paper: "Of 1,927 pQTLs, 549 (28%) were cis-acting (Supplementary Table 4). Genetic variants that change protein structure may result in apparent cis pQTLs owing to altered aptamer binding rather than true quantitative differences in protein levels. We found evidence against such artefactual associations for 371 (68%) cis pQTLs (Methods, Supplementary Tables 4, 7, 8). The results were materially unchanged when we repeated downstream analyses but excluded pQTLs without evidence against binding effects."
      - Methods section title "Evidence against aptamer-binding effects at cis pQTLs".
    - What significance threshold should we use?
    - What to do when there is cross-reactivity?
    - Some genes have multiple SOMA IDs (assayed more than once)
      - 3,195 total SOMA IDs goes to 2,872 unique Ensembl gene IDs (323 lost)

Output file naming convention:
  `{data_type}/{exp_type}/{source}/{cell_type}/{chrom}.pval{pval}.{cis_reg|trans_reg}.tsv.gz`

Output files:
  - `{chrom}.pval{pval}.cis_reg.tsv.gz`: cis-regulatory associations for `{ensembl_id}` within 1Mb (`config['sun2018_cis_window']`) of gene TSS. Filtered by pval < 2.5e-5 (`config['sun2018_cis_pval']`).
  - `{chrom}.pval{pval}.trans_reg.tsv.gz`: trans-regulatory associations for `{ensembl_id}` outside 1Mb (`config['sun2018_cis_window']`) of gene TSS. Filtered by pval < 5e-8 (`config['sun2018_trans_pval']`).

QTL output columns:
  - `chrom`, `pos`: genomic coordinates on GRCh37
  - `other_allele`: non-effect allele. Ref allele in GRCh37.
  - `effect_allele`: allele on which the effect (`beta`) is measured. Effect should always be ALT allele in GRCh37!
  - `ensembl_id`: Ensembl gene ID
  - `beta`: the effect size and direction
  - `se`: standard error of `beta`
  - `pval`: p-value of association

### Other datasets
- Closest gene
  - Use the input Ensembl VCF with VEP consequences as a variant index (`gs://genetics-portal-data/homo_sapiens_incl_consequences.vcf.gz`)
    - Convert VCF to BED using variant ID (`chrom_pos_ref_alt`) as "name"
  - Download Ensembl GRCh37 GTF gene definitions
    - Convert GTF to TSS bed file using (i) protein coding genes only, (ii) all genes.
  - Use `bedtools closest` to find closest (i) protein coding gene, (ii) any gene TSS to each variant in the index

Closest gene output columns:
  - `varid`: variant ID (`chrom_pos_ref_alt`) (GRCh37)
  - `ensemblid_protein_coding`: Ensembl ID of closest protein coding gene
  - `distance_protein_coding`: Distance (bp) from variant to closest protein coding gene
  - `ensemblid_any`: Ensembl ID of closest gene (any)
  - `distance_any`: Distance (bp) from variant to closest gene (any)

#### Requirements

Requries:
  - (Linux) split v>=8.23
  - (Mac) gsplit v>=8.23

#### Usage

```
# Install dependencies into isolated environment
conda env create -n v2g_data --file environment.yaml

# Activate environment
source activate v2g_data

# Alter configuration file
nano config.yaml

# Authenticate google cloud storage
gcloud auth application-default login

# Make manifests for QTL datasets
snakemake -s scripts/sun2018_pqtl.make_manifest.Snakefile

# Execute workflows (locally)
ncores=3
snakemake -s sun2018_pqtl.Snakefile --cores $ncores
snakemake -s gtex7_eqtl.Snakefile --cores $ncores
snakemake -s andersson2014_fantom5.Snakefile --cores $ncores
snakemake -s thurman2012_dhscor.Snakefile --cores $ncores
snakemake -s javierre2016_pchic.Snakefile --cores $ncores

# Warning. The closest gene pipeline takes ~3 hours.
snakemake -s closest_gene.Snakefile --cores $ncores --resources threads=$ncores

# Copy output to GCS
gsutil -m rsync -r -x ".*DS_Store$" output gs://genetics-portal-staging/v2g
```

#### Notes

None
