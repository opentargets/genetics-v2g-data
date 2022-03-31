Variant-to-gene functional datasets
===================================

Workflows to generate cis-regulatory datasets used for variant-to-gene (V2G) assignment in Open Targets Genetics.

## Setup

Requires:
  - Conda
  - (Linux) split v>=8.23
  - (Mac) gsplit v>=8.23

If necessary, install conda.
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
echo export PATH="$HOME/miniconda/bin:\$PATH" >> ~/.profile
. ~/.profile
```

## Generate chromatin interaction datasets
This only needs to be run if these datasets have been updated.

```bash
# Install dependencies into isolated environment
conda env create -n v2g_data --file environment.yaml

# Activate environment
conda activate v2g_data

# Alter configuration file
nano configs/config.yaml

# Authenticate google cloud storage
gcloud auth application-default login

# Increase spark memory if running locally
export PYSPARK_SUBMIT_ARGS="--driver-memory 8g pyspark-shell"

# Execute workflows (locally)
ncores=3
snakemake -s andersson2014_fantom5.Snakefile --cores $ncores
snakemake -s thurman2012_dhscor.Snakefile --cores $ncores
# snakemake -s javierre2016_pchic.Snakefile --cores $ncores
snakemake -s jung2019_pchic.Snakefile --cores $ncores

# Copy output to GCS
gsutil -m rsync -r -x ".*DS_Store$" output gs://genetics-portal-staging/v2g

```

## Process molecular QTL datasets
This needs to be run if QTL datasets have been added or updated.

```bash
# Update list of QTL studies in process_QTL_datasets_from_sumstats.w_hack.py

# Start dataproc server
gcloud beta dataproc clusters create \
    em-qtlprocess \
    --region europe-west1 \
    --image-version=preview \
    --properties=spark:spark.debug.maxToStringFields=100,spark:spark.executor.cores=32,spark:spark.executor.instances=1 \
    --master-machine-type=n2-standard-32 \
    --master-boot-disk-size=2TB \
    --num-master-local-ssds=0 \
    --zone=europe-west1-d \
    --initialization-action-timeout=20m \
    --single-node \
    --max-idle=10m

# Submit QTL processing job
gcloud dataproc jobs submit pyspark \
    --cluster=em-qtlprocess \
    --region europe-west1 \
    process_QTL_datasets_from_sumstats.w_hack.py

# To monitor
gcloud compute ssh em-qtlprocess-m \
  --project=open-targets-genetics-dev \
  --zone=europe-west1-d -- -D 1080 -N
"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
  --proxy-server="socks5://localhost:1080" \
  --user-data-dir="/tmp/em-qtlprocess-m" http://em-qtlprocess-m:8088
```

## Make mappings file

Follow steps in mapping/README.md, if new datasets have been added or updated.

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

#### eQTL

Most data now comes from eQTL catalogue, including GTEx eQTLs, which have been reprocessed by eQTL catalogue.

- [Publication link](https://pubmed.ncbi.nlm.nih.gov/34493866/)
- Dataset: summary statistics ingested from eQTL catalogue and other eQTL datasets
- Score = -log10(pval)
- Filters:
  - Variants are filtered to keep those with p-value ≤ (0.05 / number of variants tested for the gene)
  - Pre-filtered at source to remove associations >1Mb from gene TSS

#### sQTL

This is still to be done.
sQTL data will be integrated initially from GTEx v8, but may be switched to uniform data from the eQTL catalogue in the future, once this is ready.

- [Publication link](https://www.ncbi.nlm.nih.gov/pubmed/29022597)
- Score = -log10(pval)
- Filters:
  - Same as for eQTLs
#### pQTL (Sun *et al.*, 2018)

Evidence linking genetic variation to protein abundance in Sun *et al.* (2018) pQTL data.
- [Publication link](https://www.ncbi.nlm.nih.gov/pubmed/29875488)
- [Dataset link](http://www.phpc.cam.ac.uk/ceu/proteins/)
- Score = -log10(pval)
- Filters:
  - Variants are filtered to keep those with p-value ≤ (0.05 / number of variants tested for the gene)
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

#### Promoter Capture Hi-C (Jung, 2019)

Evidence linking genetic variation to genes using Promoter Capture Hi-C in 27 human cell/tissue types.
- [Publication link](https://www.nature.com/articles/s41588-019-0494-8)
- [Dataset link: Table S3](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0494-8/MediaObjects/41588_2019_494_MOESM3_ESM.xlsx)
- Score = 1 or 0
- Filters:
  - Only significant P–O pcHi-C interactions included at source

#### Promoter Capture Hi-C (Javierre, 2016)

Evidence linking genetic variation to genes using Promoter Capture Hi-C in each of the 17 human primary hematopoietic cell types.
- [Publication link](https://www.ncbi.nlm.nih.gov/pubmed/27863249)
- [Dataset link](ftp://ftp.ebi.ac.uk/pub/contrib/pchic/CHiCAGO/)
- Score = CHiCAGO score (Cairns et al., 2016)
- Filters:
  - Pre-filtered at source to contain only CHiCAGO >= 5
  - Filtered interactions across chromosomes and over 2.45e6 bases away on either side (see [plots](https://github.com/opentargets/v2g_data/blob/master/analysis/pchic%20scores%20and%20distance.ipynb))

#### Enhancer-TSS correlation (FANTOM5)

Evidence linking genetic variation to genes using correlation between the transcriptional activity of enhancers and transcription start sites using the FANTOM5 CAGE expression atlas.
- [Publication link](https://www.ncbi.nlm.nih.gov/pubmed/24670763)
- [Dataset link](http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed)
- Score = R-squared
- Filters:
  - Pre-filtered at source to contain only false-discovery rate (FDR) > 0.05
  - Pre-filtered at source to remove interactions >2Mb
  - Filtered interactions between chromosomes and over 1,068,265 bases (2 stdev) on either side (see [plots](https://github.com/opentargets/v2g_data/blob/master/analysis/fantom5.ipynb))

#### DHS-promoter correlation (Thurman, 2012)

Evidence linking genetic variation to genes using correlation of DNase I hypersensitive site and gene promoters across 125 cell and tissue types from ENCODE.
- [Publication link](https://www.ncbi.nlm.nih.gov/pubmed/22955617)
- [Dataset link](http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/dhs_gene_connectivity/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz)
- Score = R-squared
- Filters:
  - Pre-filtered at source to contain R-squared > 0.7 only
  - Pre-filtered at source to remove interactions over 500 kb
  - Pre-filtered at source to between chromsome interactions


#### (OLD) Copy from staging to genetics-portal-data

Commands for copying from staging.
NOTE: This is normally taken care of by the backend team, and the code below hasn't been updated.

```
version_date=`date +%y%m%d`

# Interval datasets
gsutil -m rsync -r gs://genetics-portal-staging/v2g/interval/pchic/javierre2016/$version_date/ gs://genetics-portal-data/v2g/$version_date/interval/pchic/javierre2016/
gsutil -m rsync -r gs://genetics-portal-staging/v2g/interval/fantom5/andersson2014/$version_date/ gs://genetics-portal-data/v2g/$version_date/interval/fantom5/andersson2014/
gsutil -m rsync -r gs://genetics-portal-staging/v2g/interval/dhscor/thurman2012/$version_date/ gs://genetics-portal-data/v2g/$version_date/interval/dhscor/thurman2012/

# QTL datasets
gsutil -m rsync -r gs://genetics-portal-staging/v2g/qtl/eqtl/gtex_v7/$version_date/ gs://genetics-portal-data/v2g/$version_date/qtl/eqtl/gtex_v7/
gsutil -m rsync -r gs://genetics-portal-staging/v2g/qtl/pqtl/sun2018/$version_date/ gs://genetics-portal-data/v2g/$version_date/qtl/pqtl/sun2018/
```
