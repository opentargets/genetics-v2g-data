#!/bin/sh

set -euo pipefail

# Args
ncores=15

# Load environment
source activate v2g_data
export PYSPARK_SUBMIT_ARGS="--driver-memory 8g pyspark-shell"

# Make manifests for QTL datasets
snakemake -s scripts/sun2018_pqtl.make_manifest.Snakefile

# Execute workflows
snakemake -s sun2018_pqtl.Snakefile --cores $ncores
snakemake -s gtex7_eqtl.Snakefile --cores $ncores
snakemake -s andersson2014_fantom5.Snakefile --cores $ncores
snakemake -s thurman2012_dhscor.Snakefile --cores $ncores
snakemake -s javierre2016_pchic.Snakefile --cores $ncores

# Copy output to GCS
gsutil -m rsync -r -x ".*DS_Store$" output gs://genetics-portal-staging/v2g

echo COMPLETE
