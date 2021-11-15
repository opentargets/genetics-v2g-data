Make biofeature mappings for QTL and other datasets
===================================================

The backend needs a table that maps from each QTL dataset to the tissue ontology code, text label, etc. Here we build that table.

Note: This has been quite confusing recently. It seems that "biofeature_code" is the field that is used in the mapping table for V2G datasets such as PC-Hi-C, so this field also needs to be what the QTL data use as the mapping for their bio_feature column, which the QTL data are partitioned on. In the past this has been a mix of either a name/label (e.g. "blood") or an ontology code (e.g. "UBERON_0000178").
But in summary, the biofeature field used in the QTL data should match the biofeature_code in the mappings file.

## Create mapping table for eQTL catalogue
```
# For eQTL catalogue, QTL metadata comes from here
wget https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv -O eqtl_catalogue.tabix_ftp_paths.tsv

# Currently we only use gene expression sumstats, not other QTL types
head -n 1 eqtl_catalogue.tabix_ftp_paths.tsv > eqtl_catalogue.tabix_ftp_paths.flt.tsv
cat eqtl_catalogue.tabix_ftp_paths.tsv | awk 'BEGIN {FS="\t"} $7 == "ge" || $7 == "microarray"' >> eqtl_catalogue.tabix_ftp_paths.flt.tsv

# Convert the table to json format with the relevant fields
python get_eqtl_catalogue_mappings.py
```

## Merge together mapping tables from separate sources

Mappings files for individual datasets, e.g. Jung 2019, Sun 2018, were made manually.

```
cat jung2019.mappings.json \
    javierre2016.mappings.json \
    pqtl_sun2018.mappings.json \
    eqtlgen_2018.mappings.json \
    eqtl_catalogue.mappings.json \
    > biofeature_labels.json
```

These mappings are used in the v2g pipeline to match QTL variants to biofeatures/tissue codes.
## Upload to GCS
```
gsutil cp biofeature_labels.json gs://genetics-portal-dev-staging/lut/211115/biofeature_labels.json
```
