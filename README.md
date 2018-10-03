Get cis-regulatory datasets
===========================

Workflows to generate cis-regulatory datasets used for variant-to-gene (V2G) assignment in Open Targets Genetics.

## Contents

- [Requirements](#requirements)
- [Usage](#usage)
- Datasets
  - [QTL](#quantitative-trait-loci-qtl-datasets)
    - [eQTL (GTEx V7)](#eqtl-gtex-v7)
    - [pQTL (Sun *et al.*, 2018)](#pqtl-sun-et-al-2018)
  - [Interval (interaction based)](#interval-interaction-based-datasets)
    - [PCHi-C (Javierre, 2016)](#promoter-capture-hi-c-javierre-2016)
    - [Enhancer-TSS corr (FANTOM5)](#enhancer-tss-corr-fantom5)
    - [DHS-promoter corr (Thurman, 2012)](#dhs-promoter-corr-thurman-2012)
  - [Functional prediction](#functional-prediction-datasets)
    - [VEP (Ensembl)](#variant-effect-predictor-vep)
  - [Closest gene](#closest-gene)

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

## Datasets

### Quantitative trait loci (QTL) datasets

QTL datatypes contain association statistics between individual genetic variation and a molecular phenotype (e.g. expression QTLs, and protein QTLs).

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

#### eQTL (GTEx V7)

Evidence linking genetic varition to gene expression in each of the 44 GTEx V7 tissues.
- [Publication link](https://www.ncbi.nlm.nih.gov/pubmed/29022597)
- Dataset: GTEx significant pairs files
- Filters:
  - Pre-filtered at source to contain only "significant" associations based on permutation analysis (p_perm < 0.05)
  - Pre-filtered at source to remove associations >1Mb from gene TSS

#### pQTL (Sun *et al.*, 2018)

Evidence linking genetic varition to protein abundance in Sun *et al.* (2018) pQTL data.
- [Publication link](https://www.ncbi.nlm.nih.gov/pubmed/29875488)
- [Dataset link](http://www.phpc.cam.ac.uk/ceu/proteins/)
- Filters:
  - Filtered to approximately match GTEx permutation filtering using a p < 2.5e-5 threshold. This equates to ~2000 independant tests per gene region.
  - Filtered to remove associations >1Mb from each gene TSS.
- Notes:
  - 3,195 total SOMA IDs goes to 2,872 unique Ensembl gene IDs (323 lost). For these we are deduplicating by gene_id and keeping the first occurence in the manifest.
  - We do not exclude "pQTLs without evidence against binding effects".

 > *"Of 1,927 pQTLs, 549 (28%) were cis-acting (Supplementary Table 4). Genetic variants that change protein structure may result in apparent cis pQTLs owing to altered aptamer binding rather than true quantitative differences in protein levels. We found evidence against such artefactual associations for 371 (68%) cis pQTLs (Methods, Supplementary Tables 4, 7, 8). The results were materially unchanged when we repeated downstream analyses but excluded pQTLs without evidence against binding effects."* Sun *et al.*

### Interval (interaction based) datasets

Interval datatypes contain genomic regions that are linked to gene transcription start site through molecular interactions (e.g. promoter Capture Hi-C).

Output file naming convention:
  `{data type}/{experiment type}/{source}/{cell|tissue type}/{chrom}.{processed|raw}.tsv.gz`

Interval output columns:
  - `chrom`, `start`, `end`: 1-based genomic co-ordinates. Start and end are inclusive.
  - `ensembl_id`: Ensembl gene ID
  - `score`: feature score from original dataset (unmodified)
  - `cell_type`: cell or tissue type of this experiment

#### Promoter Capture Hi-C (Javierre, 2016)

Evidence linking genetic variation to genes using Promoter Capture Hi-C in each of the 17 human primary hematopoietic cell types.
- [Publication link](https://www.ncbi.nlm.nih.gov/pubmed/27863249)
- [Dataset link](ftp://ftp.ebi.ac.uk/pub/contrib/pchic/CHiCAGO/)
- Score = CHiCAGO score (Cairns et al., 2016)
- Filters:
  - Pre-filtered at source to contain only CHiCAGO >= 5
  - Filtered interactions across chromosomes and over 2.45e6 bases away on either side (see [plots](https://github.com/opentargets/v2g_data/blob/master/analysis/pchic%20scores%20and%20distance.ipynb))

#### Enhancer-TSS corr (FANTOM5)


Evidence linking genetic variation to genes using correlation between the transcriptional activity of enhancers and transcription start sites using the FANTOM5 CAGE expression atlas.
- [Publication link](https://www.ncbi.nlm.nih.gov/pubmed/24670763)
- [Dataset link](http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed)
- Score = R-squared
- Filters:
  - Pre-filtered at source to contain only false-discovery rate (FDR) > 0.05
  - Pre-filtered at source to remove interactions >2Mb
  - Filtered interactions between chromosomes and over 1,068,265 bases (2 stdev) on either side (see [plots](https://github.com/opentargets/v2g_data/blob/master/analysis/fantom5.ipynb))

#### DHS-promoter corr (Thurman, 2012)

Evidence linking genetic variation to genes using correlation of DNase I hypersensitive site and gene promoters across 125 cell and tissue types from ENCODE.
- [Publication link](https://www.ncbi.nlm.nih.gov/pubmed/22955617)
- [Dataset link](http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/dhs_gene_connectivity/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz)
- Score = R-squared (R-squared > 0.7 only)
- Filters:
  - Pre-filtered at source to contain R-squared > 0.7 only
  - Pre-filtered at source to remove interactions over 500 kb
  - Pre-filtered at source to between chromsome interactions

### Functional prediction datasets

Functional prediction datatypes contain *in silico* prediction of the functional consequence of a genetic variant on a gene's transcript (e.g. VEP).

#### Variant effect predictor (VEP)

Most severe coding-sequence consequence(s) from Ensembl's Variant Effect Predictor.
- [Publication link](https://www.ncbi.nlm.nih.gov/pubmed/27268795)
- [Dataset link](ftp://ftp.ensembl.org/pub/grch37/update/variation/vcf/homo_sapiens/Homo_sapiens_incl_consequences.vcf.gz)
- Scores = consequences are mapped to `v2g_score` in [this file](configs/vep_consequences.tsv)
- Filters:
  - Consequences with score 0 are removed from V2G assignment

### Closest gene

Closest gene files contain the closest gene for each variant in the variant index.

Closest gene output columns:
  - `varid`: variant ID (`chrom_pos_ref_alt`) (GRCh37)
  - `ensemblid_protein_coding`: Ensembl ID of nearest protein coding gene
  - `distance_protein_coding`: Distance (bp) from variant to nearest protein coding gene
  - `ensemblid_any`: Ensembl ID of nearest gene (any)
  - `distance_any`: Distance (bp) from variant to nearest gene (any)

Closest gene methods
  - Use the input Ensembl VCF with VEP consequences as a variant index (`gs://genetics-portal-data/homo_sapiens_incl_consequences.vcf.gz`)
    - Convert VCF to BED using variant ID (`chrom_pos_ref_alt`) as "name"
  - Download Ensembl GRCh37 GTF gene definitions
    - Convert GTF to TSS bed file using (i) protein coding genes only, (ii) all genes.
  - Use `bedtools closest` to find nearest (i) protein coding gene, (ii) any gene TSS to each variant in the index


#### Copy from staging to genetics-portal-data

Dry-run commands for copying from staging

```
version_date=`date +%y%m%d`

# nearest gene
gsutil -m rsync -rn gs://genetics-portal-staging/v2g/nearest_gene/$version_date/ gs://genetics-portal-data/v2g/nearest_gene/

# Interval datasets
gsutil -m rsync -rn gs://genetics-portal-staging/v2g/interval/pchic/javierre2016/$version_date/ gs://genetics-portal-data/v2g/interval/pchic/javierre2016/
gsutil -m rsync -rn gs://genetics-portal-staging/v2g/interval/fantom5/andersson2014/$version_date/ gs://genetics-portal-data/v2g/interval/fantom5/andersson2014/
gsutil -m rsync -rn gs://genetics-portal-staging/v2g/interval/dhscor/thurman2012/$version_date/ gs://genetics-portal-data/v2g/interval/dhscor/thurman2012/

# QTL datasets
gsutil -m rsync -rn gs://genetics-portal-staging/v2g/qtl/eqtl/gtex_v7/$version_date/ gs://genetics-portal-data/v2g/qtl/eqtl/gtex_v7/
gsutil -m rsync -rn gs://genetics-portal-staging/v2g/qtl/pqtl/sun2018/$version_date/ gs://genetics-portal-data/v2g/qtl/pqtl/sun2018/
```
