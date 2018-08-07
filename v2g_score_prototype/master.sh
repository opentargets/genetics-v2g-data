#!/usr/bin/env bash
#

#
# Scoring using V2G data -------------------------------------------------------
#

# Input v2g data
in_slices=input_v2g_slices/*.tsv.gz
# in_slices=input_v2g_slices/*FTO*.tsv.gz

# Make output dir
mkdir -p output_v2g

# Run scoring on v2g data
for slice in $in_slices; do
  # Make output name
  bn=$(basename $slice)
  outpref=output_v2g/${bn/.neighboorhood.tsv.gz/}
  # Get varid from basename
  varid=$(echo $bn | cut -f 3 -d ".")
  # Make command
  echo python scripts/score_v2g_v2.py \
    --v2g_slice $slice \
    --varid $varid \
    --outpref $outpref
done | parallel -j 3

#
# Scoring using D2V2G data -------------------------------------------------------
#

# TODO
