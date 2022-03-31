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

Mappings files for individual datasets, e.g. Jung 2019, pqtls, were made manually.

```
cd mapping
cat jung2019.mappings.json \
    javierre2016.mappings.json \
    pqtl_mappings.json \
    eqtlgen_2018.mappings.json \
    eqtl_catalogue.mappings.json \
    gtex_sqtl.mappings.json \
    > biofeature_labels.json
```

## Make composite biofeature mappings

This is a "hack" that is used by the frontend UI to display both study ID and biofeatures, grouping by biofeature.
```
python get_composite_mappings.py --in biofeature_labels.json --out biofeature_labels.composites.json
```

## Merge in composite mappings

```
cat biofeature_labels.json \
    biofeature_labels.composites.json \
    > biofeature_labels.w_composites_new.json
```

Edit this file as needed to make labels more readable, etc.
This is not an ideal process, because when you re-run the above then you will overwrite the manual edits.
What I do is to manually take only the new lines from the above process.

## Upload to GCS
```
gsutil cp biofeature_labels.w_composites.json gs://genetics-portal-dev-staging/lut/biofeature_labels/220212/biofeature_labels.w_composites.json
```
