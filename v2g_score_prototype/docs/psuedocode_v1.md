Our cis-regulatory data is structured as:

`gs://genetics-portal-staging/v2g/{data type}/{experiment type}/{source}/{version}/{tissue or cell line}/{chroms}.processed.tsv.gz`

E.g. here's one tissue from each of our data sources
```
gs://genetics-portal-staging/v2g/interval/dhscor/thurman2012/180703/unspecified/1-23.processed.tsv.gz
gs://genetics-portal-staging/v2g/interval/fantom5/andersson2014/180703/unspecified/1-23.processed.tsv.gz
gs://genetics-portal-staging/v2g/interval/pchic/javierre2016/180703/Endothelial_precursors/1-23.processed.tsv.gz
gs://genetics-portal-staging/v2g/qtl/eqtl/gtex_v7/180703/Adipose_Subcutaneous/1-23.cis_reg.processed.tsv.gz
gs://genetics-portal-staging/v2g/qtl/pqtl/sun2018/180703/Blood_plasma/1-22.pval2.5e-05.cis_reg.processed.tsv.gz
```

We would like to maintain our ability to quickly add new datasets without changing the scoring algorithm too much (except maybe adding new weights). We therefore want something generalisable and non-parametric, as scores between the different experiment types are not necessarily comparable.

```
# Pseudocode

This is assuming we've fixed for a single variant. In practice, we'll calculate for all variants in a region.

## Preprocessing (applied across whole dataset)
- Create a score column
  * QTL datasets: 1 - pvalue
  * Interval datasets: interval score
  * VEP: map to https://github.com/opentargets/postgap-api/blob/master/data-prep/eco_scores.tsv
- Group by (source, tissue), transform into quantiles

## Aggregate across tissues (per source)
- Group by (source, gene)
- Calculate features: (i) max score, (ii) non-null tissue count
- Scale feature (ii) to 0-1, so that its on a similar scale to (i)
- Take weighted mean of features (i) and (ii) using ft_weights (see below)
- This gives a score per gene for each source

## Aggregate across sources (per variant)
- Group by (gene)
- Caluclate weighted mean over sources using source weights (see below)
- This gives a score per gene for a given variant

ft_weights = {
    'vep':           {'ft_max': 1,   'ft_count_scaled': 0},
    'javierre2016':  {'ft_max': 0.8, 'ft_count_scaled': 0.2},
    'andersson2014': {'ft_max': 0.8, 'ft_count_scaled': 0.2},
    'thurman2012':   {'ft_max': 0.8, 'ft_count_scaled': 0.2},
    'gtex_v7':       {'ft_max': 0.8, 'ft_count_scaled': 0.2},
    'sun2018':       {'ft_max': 0.8, 'ft_count_scaled': 0.2}
}

source_weights = {
    'vep': 1,
    'javierre2016': 1,
    'andersson2014': 1,
    'thurman2012': 1,
    'gtex_v7': 1,
    'sun2018': 1
}

```

As mentioned above, we may actually want to calculate a 2nd score (locus-to-gene) which will be used for the main platform (as well as places in the genetics portal). In this, for a given index variant we will aggregate over tag varinats in the LD or finemapping set using another weighted mean.
