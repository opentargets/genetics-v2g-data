#!/usr/bin/env bash
#

set -euo pipefail

clickhouse-client -q "select distinct(gene_id) from ot.d2v2g" > d2v2g_gene_list.tsv
clickhouse-client -q "select distinct(gene_id) from ot.v2g" > v2g_gene_list.tsv

echo COMPLETE
