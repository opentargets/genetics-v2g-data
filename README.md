Get cis-regulatory datasets
===========================

Process cis-regulatory datasets for variant-to-gene (V2G) assignment.

#### Datasets

### Interval datasets
- PCHiC from Javierre et al 2016
  - [Publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5123897/)
  - [Dataset](ftp://ftp.ebi.ac.uk/pub/contrib/pchic/CHiCAGO/))
  - Score = CHiCAGO score (Cairns et al., 2016) (CHiCAGO >= 5 only)
  - Filtering interactions across chromosomes and over 2.45e6 bases away on either side (see [plots](https://github.com/opentargets/v2g_data/blob/master/analysis/pchic%20scores%20and%20distance.ipynb))
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
  - cis-regulatory "significant" associations only
- pQTL from Sun et al 2018
  - [Publication](https://www.nature.com/articles/s41586-018-0175-2)
  - [Data available from](http://www.phpc.cam.ac.uk/ceu/proteins/)
  - Data downloaded to `gs://genetics-portal-raw/pqtl_sun2018`
  - Maunipulated to be similar to GTEx datasets:
    - Keeping only cis-regulatory (within 1 Mb of gene TSS)
    - filtered to p < 2.5e-5 to approximately match GTEx
  - Notes:
    - Do we only keep cis-associations?
    - Do we need to exclude "pQTLs without evidence against binding effects"
      - Main body of paper: "Of 1,927 pQTLs, 549 (28%) were cis-acting (Supplementary Table 4). Genetic variants that change protein structure may result in apparent cis pQTLs owing to altered aptamer binding rather than true quantitative differences in protein levels. We found evidence against such artefactual associations for 371 (68%) cis pQTLs (Methods, Supplementary Tables 4, 7, 8). The results were materially unchanged when we repeated downstream analyses but excluded pQTLs without evidence against binding effects."
      - Methods section title "Evidence against aptamer-binding effects at cis pQTLs".
    - What to do when there is cross-reactivity?
    - Some genes have multiple SOMA IDs (assayed more than once)
      - 3,195 total SOMA IDs goes to 2,872 unique Ensembl gene IDs (323 lost)
      - For these we are deduplicating by gene_id and keeping the first occurence

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
- Nearest gene
  - Use the input Ensembl VCF with VEP consequences as a variant index (`gs://genetics-portal-data/homo_sapiens_incl_consequences.vcf.gz`)
    - Convert VCF to BED using variant ID (`chrom_pos_ref_alt`) as "name"
  - Download Ensembl GRCh37 GTF gene definitions
    - Convert GTF to TSS bed file using (i) protein coding genes only, (ii) all genes.
  - Use `bedtools closest` to find nearest (i) protein coding gene, (ii) any gene TSS to each variant in the index

Nearest gene output columns:
  - `varid`: variant ID (`chrom_pos_ref_alt`) (GRCh37)
  - `ensemblid_protein_coding`: Ensembl ID of nearest protein coding gene
  - `distance_protein_coding`: Distance (bp) from variant to nearest protein coding gene
  - `ensemblid_any`: Ensembl ID of nearest gene (any)
  - `distance_any`: Distance (bp) from variant to nearest gene (any)

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

# Warning. The nearest gene pipeline takes ~3 hours.
snakemake -s nearest_gene.Snakefile --cores $ncores --resources threads=$ncores

# Copy output to GCS
gsutil -m rsync -r -x ".*DS_Store$" output gs://genetics-portal-staging/v2g
```

#### Copy from staging to genetics-portal-data

Dry-run commands for copying from staging

```
# nearest gene
gsutil -m rsync -rn gs://genetics-portal-staging/v2g/nearest_gene/180830/ gs://genetics-portal-data/v2g/nearest_gene/

# Interval datasets
gsutil -m rsync -rn gs://genetics-portal-staging/v2g/interval/pchic/javierre2016/180831/ gs://genetics-portal-data/v2g/interval/pchic/javierre2016/
gsutil -m rsync -rn gs://genetics-portal-staging/v2g/interval/fantom5/andersson2014/180831/ gs://genetics-portal-data/v2g/interval/fantom5/andersson2014/
gsutil -m rsync -rn gs://genetics-portal-staging/v2g/interval/dhscor/thurman2012/180831/ gs://genetics-portal-data/v2g/interval/dhscor/thurman2012/

# QTL datasets
gsutil -m rsync -rn gs://genetics-portal-staging/v2g/qtl/eqtl/gtex_v7/180809/ gs://genetics-portal-data/v2g/qtl/eqtl/gtex_v7/
gsutil -m rsync -rn gs://genetics-portal-staging/v2g/qtl/pqtl/sun2018/180808/ gs://genetics-portal-data/v2g/qtl/pqtl/sun2018/
```
