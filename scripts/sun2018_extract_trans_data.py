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

    # Load each chromosome in chunks to save memory
    dfs = []
    for inf in args.infs:

        iter_tsv = pd.read_csv(inf, sep='\t', header=0, iterator=True,
                               chunksize=100000)
        for chunk in iter_tsv:

            # Format data types
            chunk.chromosome = chunk.chromosome.astype(str)
            chunk.position = chunk.position.astype(int)
            chunk.Allele1 = chunk.Allele1.str.upper()
            chunk.Allele2 = chunk.Allele2.str.upper()

            # Filter by p
            chunk['pval'] = chunk['log(P)'].rpow(10)
            chunk = chunk.loc[chunk.pval <= args.pval, :]

            # If on the same chromosome as gene, select variants outside the
            # cis-regulatory region. If on different chromosome, select everything.
            to_keep = (((chunk['chromosome'] == args.chrom) &
                        ((chunk['position'] <= (args.tss_pos - args.cis_window)) |
                        (chunk['position'] >= (args.tss_pos + args.cis_window))))
                       | (chunk['chromosome'] != args.chrom))
            chunk = chunk.loc[to_keep, :]

            dfs.append(chunk)

    # Merge chunks
    data = pd.concat(dfs)

    # Select required columns
    data = data.rename(columns={
        'chromosome': 'chrom',
        'position': 'pos',
        'Allele1': 'effect_allele',
        'Allele2': 'other_allele',
        'Effect': 'beta',
        'StdErr': 'se',
        'pval': 'pval'})

    data = data.loc[:, ['chrom','pos', 'other_allele', 'effect_allele',
                        'beta', 'se', 'pval']]
    data = data.sort_values(['chrom', 'pos'])

    data.to_csv(args.outf, sep='\t', index=None, compression='gzip')

    return 0

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--infs', metavar="<file>", help=('Input file'), type=str, nargs='+', required=True)
    parser.add_argument('--outf', metavar="<file>", help=("Output file"), type=str, required=True)
    parser.add_argument('--chrom', metavar="<str>", help=("Chromosome"), type=str, required=True)
    parser.add_argument('--tss_pos', metavar="<int>", help=("TSS position"), type=int, required=True)
    parser.add_argument('--cis_window', metavar="<int>", help=("window around TSS to extract data from"), type=int, required=True)
    parser.add_argument('--pval', metavar="<float>", help=("Minimum P-value to be considered significant"), type=float, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
