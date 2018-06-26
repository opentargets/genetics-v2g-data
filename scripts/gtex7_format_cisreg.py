#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Ed Mountjoy (June 2018)

Formats the output of bedtools into final format for loading
'''

import argparse
import pandas as pd
import sys


def main():

    # Parse args
    args = parse_args()

    # Load data
    data = pd.read_csv(args.inf, sep='\t', header=0)

    # Split variant ID into parts
    var_split = data['variant_id'].str.split('_', expand=True)
    var_split.columns = ['chrom', 'pos', 'other_allele', 'effect_allele', 'build']

    # Join
    data = pd.concat([var_split, data], axis=1)
    data['chrom'] = data['chrom'].astype(str)
    data['pos'] = data['pos'].astype(int)

    # Fix gene_ids
    data['gene_id'] = data['gene_id'].apply(lambda gene: gene.split('.')[0])

    # Replace chrom X for 23
    # data = data.replace(to_replace={'chrom':{'X': '23'}})

    # Select required columns
    data = data.rename(columns={
        'slope': 'beta',
        'slope_se': 'se',
        'gene_id': 'ensembl_id',
        'pval_nominal': 'pval'})
    data = data.loc[:, ['chrom','pos', 'other_allele', 'effect_allele',
                        'ensembl_id', 'beta', 'se', 'pval']]
    data = data.sort_values(['chrom', 'pos'])

    data.to_csv(args.outf, sep='\t', index=None, compression='gzip')

    return 0

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('Input file'), type=str, required=True)
    parser.add_argument('--outf', metavar="<file>", help=("Output file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
