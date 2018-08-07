```
# Pseudocode

This is assuming we've fixed for a single variant. In practice, we'll calculate for all variants in a region.

## Preprocessing (applied across whole dataset)
- Create a score column
  * QTL datasets: 1 - pvalue
  * Interval datasets: interval score
  * VEP map to v2g_score column in gs://genetics-portal-data/lut/vep_consequences.tsv
- Group by (source, tissue), transform into quantiles (I suggest 10 bins [deciles] but this is flexible)

## Aggregate across tissues (per source)
- Group by (source, gene)
- Calculate max score across tissues
- This gives a score per gene for each source

## Aggregate across sources (per variant)
- Group by (gene)
- Caluclate weighted mean over sources using source weights (see below)
- This gives a score per gene for a given variant

source_weights = {
    'vep': 1,
    'javierre2016': 1,
    'andersson2014': 1,
    'thurman2012': 1,
    'gtex_v7': 1,
    'sun2018': 1
}

```
