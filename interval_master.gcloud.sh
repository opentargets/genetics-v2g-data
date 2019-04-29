#!/bin/sh

set -euo pipefail

# Args
ncores=15
instance_zone=europe-west1-d
instance_name=em-v2g

# Load environment
source activate v2g_data
export PYSPARK_SUBMIT_ARGS="--driver-memory 8g pyspark-shell"

# Execute workflows
snakemake -s andersson2014_fantom5.Snakefile --cores $ncores
snakemake -s thurman2012_dhscor.Snakefile --cores $ncores
snakemake -s javierre2016_pchic.Snakefile --cores $ncores

# Copy output to GCS
gsutil -m rsync -r -x ".*DS_Store$" output gs://genetics-portal-staging/v2g
gcloud compute instances stop $instance_name --zone=$instance_zone

echo COMPLETE
