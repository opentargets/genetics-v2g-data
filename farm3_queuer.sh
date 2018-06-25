#!/bin/sh
#BSUB -J get_g2v_data
#BSUB -q long
#BSUB -n 8
#BSUB -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" -M16000
#BSUB -o output.%J
#BSUB -e errorfile.%J

source activate g2v_data
export GCLOUD_PROJECT=open-targets-genetics

snakemake --cores 7
