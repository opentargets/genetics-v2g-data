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

    # Format data types
    data.chromosome = data.chromosome.astype(str)
    data.position = data.position.astype(int)
    data.Allele1 = data.Allele1.str.upper()
    data.Allele2 = data.Allele2.str.upper()

    # Convert to p
    data['pval'] = data['log(P)'].rpow(10)

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
    parser.add_argument('--inf', metavar="<file>", help=('Input file'), type=str, required=True)
    parser.add_argument('--outf', metavar="<file>", help=("Output file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
