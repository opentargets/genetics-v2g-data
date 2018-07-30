#!/usr/bin/env bash
#

in_slices=data/*.csv
in_slices=data/*.FTO.csv

for slice in $in_slices; do
  bn=$(basename $slice)
  python scripts/score_test_v1.py \
    --v2g_slice $slice \
    --out output/${bn/.neighboorhood.csv/.scores.txt}
done
