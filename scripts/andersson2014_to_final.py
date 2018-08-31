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
    # Map gene symbols to ensembl gene ids
    #

    # Make and apply map
    gene_map = make_symbol_to_ensembl_map(args.gtf)
    gene_ids = df3['gene_symbol'].apply(lambda x: gene_map.get(x.upper(), None))

    # Swap symbol for ens id
    df3.insert(3, 'ensembl_id', gene_ids)
    df3 = df3.drop('gene_symbol', axis=1)

    # Remove rows with None
    df3 = df3.loc[~pd.isnull(df3['ensembl_id']), :]

    #
    # Format output
    #

    # Sort
    df3.chrom = df3.chrom.astype(str)
    df3 = df3.sort_values(['chrom', 'start', 'end', 'score'])

    # Remove duplicates
    df3 = df3.drop_duplicates(
        subset=['chrom', 'start', 'end', 'ensembl_id'],
        keep='last')

    # Add 1 to make it 1-based coords
    df3.loc[:, 'start'] = df3['start'].astype(int)
    df3.loc[:, 'end'] = df3['end'].astype(int)
    df3['start'] = df3['start'] + 1

    # Outpt
    df3.loc[:, 'cell_type'] = args.cell_name
    df3.to_csv(args.outf, sep='\t', index=None, compression='gzip')

    return 0

def make_symbol_to_ensembl_map(gtf_file):
    ''' Parses an Ensembl GTF to extract symbol to gene id map
    Args:
        gtf_file (str): Ensembl gtf file
    Returns:
        Dictionary {symbol: id}
    '''
    d = {}
    with gzip.open(gtf_file, 'r') as in_h:
        for line in in_h:
            # Skip comments
            line = line.decode()
            if line.startswith('#'):
                continue
            # Only keep gene features
            parts = line.rstrip().split('\t')
            if not parts[2] == 'gene':
                continue
            # Extract symbol and gene id
            info_dict = parse_gtf_info_field(parts[8])
            d[info_dict['gene_name'].upper()] = info_dict['gene_id']
    return d

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
    parser.add_argument('--gtf', metavar="<str>", help=("Ensembl gtf file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
