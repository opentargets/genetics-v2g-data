#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Ed Mountjoy (June 2018)

Formats the Andersson et al tss-enhancer file to final output
'''

import sys
import argparse
import pandas as pd
import gzip

def main():

    # Parse args
    args = parse_args()

    #
    # Load and extract required columns
    #

    # Load full file
    df = pd.read_csv(args.inf, sep='\t', header=0, low_memory=False, skiprows=1)

    # Split "name" column
    df2 = df['name'].str.split(';', expand=True)
    df2 = df2.iloc[:, [0, 2, 3]]

    # Split interval column
    loc = df2.iloc[:, 0].str.split(':|-', expand=True)

    # Combine
    df3 = pd.concat([loc,
                     df2.loc[:, 2],
                     df2.loc[:, 3].apply(lambda x: float(x.replace('R:', '')))
                     ], axis=1)
    df3.columns = ["chrom", "start", "end", "gene_symbol", "score"]

    # Sqaure the score
    df3['score'] = df3['score']**2

    # Strip "chr" from chrom
    df3.loc[:, 'chrom'] = df3.loc[:, 'chrom'].str.replace('chr', '')

    #
    # Map gene symbols to gene information (ensembl ID, chromosme, tss position)
    #

    # Load gene info json
    gene_info = pd.read_json(args.gene_info, lines=True)
    gene_info = gene_info.loc[:, ['gene_name', 'gene_id', 'chr', 'tss']]

    # Merge gene info to data
    df3.gene_symbol = df3.gene_symbol.str.upper()
    df3 = pd.merge(df3, gene_info, left_on='gene_symbol', right_on='gene_name',
                   how='left')

    # Remove rows with None
    df3 = df3.loc[~pd.isnull(df3['gene_id']), :]

    #
    # Filter trans and non-proximal (>2*5.341328e05) interactions
    #

    total_lines = df3.shape[0]

    # Remove interactions between different chromosomes aka. "trans"
    cis = (df3['chrom'] == df3['chr'])
    print("CIS interactions: \t\t{} ({:.2f} %)".format(sum(cis),(sum(cis)*100)/total_lines))

    # Remove interactions where the midpoint of the interval is too far away
    # from the gene TSS
    df3['start'] = df3['start'].astype(int)
    df3['end'] = df3['end'].astype(int)
    df3['tss'] = df3['tss'].astype(int)
    df3['distance'] = abs((df3['start'] + df3['end'])/2 - df3['tss'])
    twosd_thres = 2 * 5.341328e5
    near = (df3['distance'] < twosd_thres)
    print("Interactions < {:.2E} bases: \t{} ({:.2f} %)".format(twosd_thres,sum(near),(sum(near)*100)/total_lines))

    # Apply filters and extract required columns
    filtered_df = df3.loc[cis & near]
    out_df = filtered_df.loc[:, ['chrom', 'start', 'end', 'gene_id', 'score']].copy()
    out_df.columns = ['chrom', 'start', 'end', 'ensembl_id', 'score']

    #
    # Format output
    #

    # Sort
    out_df.chrom = out_df.chrom.astype(str)
    out_df = out_df.sort_values(['chrom', 'start', 'end', 'score'])

    # Remove duplicates
    out_df = out_df.drop_duplicates(
        subset=['chrom', 'start', 'end', 'ensembl_id'],
        keep='last')

    # Add 1 to make it 1-based coords
    out_df['start'] = out_df['start'] + 1

    # Outpt
    out_df.loc[:, 'cell_type'] = args.cell_name
    out_df.to_csv(args.outf, sep='\t', index=None, compression='gzip')

    return 0

def parse_gtf_info_field(info_str):
    """ Parse gtf info string into a dictionary
    Args:
        info_str (str): info field from gtf
    Return:
        {key, value for field in info}
    """
    d = {}
    for pair in info_str.split('; '):
        key, value = pair.split(' ')
        d[key] = value.strip('"')
    return d

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('Input file'), type=str, required=True)
    parser.add_argument('--outf', metavar="<file>", help=("Output file"), type=str, required=True)
    parser.add_argument('--cell_name', metavar="<str>", help=("Name of cell type"), type=str, required=True)
    parser.add_argument('--gene_info', metavar="<str>", help=("Gene information lut"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
