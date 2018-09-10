#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Ed Mountjoy (June 2018)

Formats the output of bedtools into final format for loading
'''

import argparse
import pandas as pd

def main():

    # Parse args
    args = parse_args()

    # Load
    df = pd.read_csv(args.inf, sep='\t', header=None, low_memory=False)

    total_lines = df.shape[0]
    # Remove interactions between different chromosomes aka. "trans"
    cis = (df.iloc[:,3] == df.iloc[:,9])
    print("CIS interactions: \t\t{} ({:.2f} %)".format(sum(cis),(sum(cis)*100)/total_lines))

    # Remove interactions where the midpoint of the interval is too far away
    # from the gene TSS
    df['distance'] = abs((df.iloc[:,4] + df.iloc[:,5])/2 - df.iloc[:,10])
    twosd_thres = 2.45e6
    near = (df['distance'] < twosd_thres)
    print("Interactions < {:.2E} bases: \t{} ({:.2f} %)".format(twosd_thres,sum(near),(sum(near)*100)/total_lines))

    # Apply filters and extract required columns
    filtered_df = df.loc[cis & near]
    out_df = filtered_df.iloc[:, [3, 4, 5, 12, 6]].copy()
    out_df.columns = ['chrom', 'start', 'end', 'ensembl_id', 'score']

    # Add 1 to the start to make everything 1-indexed
    out_df.start = out_df.start + 1

    # Sort
    out_df.chrom = out_df.chrom.astype(str)
    out_df = out_df.sort_values(['chrom', 'start', 'end', 'score'])

    # Remove duplicates
    out_df = out_df.drop_duplicates(
        subset=['chrom', 'start', 'end', 'ensembl_id'],
        keep='last')

    # Write tsv
    out_df.loc[:, 'cell_type'] = args.cell_name
    out_df.to_csv(args.outf, sep='\t', index=None, compression='gzip')

    return 0

def parse_args():
    ''' Load command line args '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar='<file>', help=('Input file'), type=str, required=True)
    parser.add_argument('--outf', metavar='<file>', help=('Output file'), type=str, required=True)
    parser.add_argument('--cell_name', metavar='<str>', help=('Name of cell type'), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
