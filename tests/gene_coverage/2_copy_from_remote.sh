#!/usr/bin/env bash
#

set -euo pipefail

gcloud compute scp "emountjoy_statgen@clickhouse-0230-dev:~/*_gene_list.tsv" output

echo COMPLETE
