#!/bin/sh
#BSUB -J vep-nearest
#BSUB -q normal
#BSUB -n 1
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" -M8000
#BSUB -o output.%J
#BSUB -e errorfile.%J

# bsub -q normal -J interactive -n 1 -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" -M8000 -Is bash

# Args
vep_dir=/lustre/scratch115/projects/otcoregen/reference_data
in_gs=gs://genetics-portal-data/homo_sapiens_incl_consequences.vcf.gz

# Make input vcf
mkdir -p input_data
in_vcf=input_data/input.vcf.gz
gsutil cat $in_gs | zgrep -v "^##" | head -10000 | cut -f 1-7 | gzip -c > $in_vcf

# Run VEP
mkdir -p output
outf=output/vcf_nearest_gene.txt
$vep_dir/ensembl-vep/vep \
  --dir_cache $vep_dir/ensembl-vep \
  --cache_version 91 \
  --cache \
  --assembly GRCh37 \
  --port 3337 \
  --input_file $in_vcf \
  --format vcf \
  --output_file $outf \
  --nearest gene \
  --force_overwrite
