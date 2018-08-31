#!/bin/sh
#BSUB -J v2g_master
#BSUB -q long
#BSUB -n 8
#BSUB -R "select[mem>64000] rusage[mem=64000] span[hosts=1]" -M64000
#BSUB -o output.%J
#BSUB -e errorfile.%J

set -euo pipefail

# Args
ncores=8

# Load environment
source activate v2g_data
module load hgi/coreutils/8.23

# Make manifests for QTL datasets
snakemake -s scripts/sun2018_pqtl.make_manifest.Snakefile

# Execute workflows
snakemake -s sun2018_pqtl.Snakefile --cores $ncores
snakemake -s gtex7_eqtl.Snakefile --cores $ncores
snakemake -s andersson2014_fantom5.Snakefile --cores $ncores
snakemake -s thurman2012_dhscor.Snakefile --cores $ncores
snakemake -s javierre2016_pchic.Snakefile --cores $ncores
snakemake -s nearest_gene.Snakefile --cores $ncores --resources threads=$ncores

# Copy output to GCS
gsutil -m rsync -r -x ".*DS_Store$" output gs://genetics-portal-staging/v2g

echo COMPLETE
