#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import pandas as pd

def main():

    # Args
    in_gene_lut = '/Users/em21/Projects/ot_genetics/genetics-backend/makeLUTs/output/gene_dictionary.json'
    gene_lists = {
        'v2g_gene': 'output/v2g_gene_list.tsv',
        'd2v2g_gene': 'output/d2v2g_gene_list.tsv'
    }
    out_file = 'output/gene_v2g_membership.tsv.gz'

    # Load gene dict
    df = pd.read_json(in_gene_lut, lines=True)
    df = df.loc[df.biotype == 'protein_coding', :]

    # Add column showing if genes are in v2g
    for name, infile in gene_lists.items():

        # Load set
        gene_set = set(pd.read_csv(infile, sep='\t').iloc[:, 0])

        # Check membership
        df[name] = df.gene_id.isin(gene_set)

    # Print stats
    print('Total protein-coding genes: {0}'.format(df.shape[0]))
    print('Total protein-coding in v2g table: {0}'.format(df.v2g_gene.sum()))
    print('Total protein-coding in d2v2g table: {0}'.format(df.d2v2g_gene.sum()))

    # Save
    df.loc[:, ['biotype', 'gene_id', 'gene_name', 'description', 'v2g_gene', 'd2v2g_gene']].to_csv(out_file, sep='\t', index=None, compression='gzip')





    return 0

if __name__ == '__main__':

    main()
