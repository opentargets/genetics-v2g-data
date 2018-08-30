#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Ed Mountjoy (June 2018)

Converts Ensembl GTF to bed of Gene TSS.

Notes:
    - Keeps protein coding transcripts only.

'''

import argparse
import pandas as pd

def main():

    # Parse args
    args = parse_args()

    # Read gtf in chunks to save RAM
    iter_csv = pd.read_csv(args.inf, sep='\t', header=None, comment='#',
                           iterator=True, chunksize=100000, low_memory=False)

    # Only keep transcripts
    df = pd.concat([chunk[chunk.iloc[:, 2] == 'transcript'] for chunk in iter_csv])

    # Only keep protein coding transcripts
    if args.protein_coding:
        to_keep = df.iloc[:, 8].str.contains('transcript_biotype "protein_coding"')
        df = df.loc[to_keep, :]

    # Convert rows to bed format
    bed = df.apply(gtf_row_to_bed, axis=1)
    bed.columns = ['chrom', 'start', 'end', 'name', 'score']
    bed.chrom = bed.chrom.astype(str)
    bed.start = bed.start.astype(int)
    bed.end = bed.end.astype(int)
    bed = bed.sort_values(['chrom', 'start', 'end'])

    # Drop duplicate rows (multiple transcripts)
    bed = bed.drop_duplicates()

    bed.to_csv(args.outf, sep='\t', header=None, index=None, compression='gzip')

    return 0

def gtf_row_to_bed(row):
    """ Converts gtf row to bed format:
    Args:
        row (pd.Series): One GTF row
    Returns:
        pd.Series([chrom, start, end, ensID, "."])
    """
    # Extract required fields
    chrom = str(row[0])
    strand = row[6]
    # Take end of transcript for forward/reverse strand
    if strand == "+":
        end = int(row[3])
    elif strand == "-":
        end = int(row[4])
    # Bed files are 0-indexed (whereas GTF are 1-index) so start should be end-1
    start = end - 1
    ensID = parse_gtf_info_field(row[8])['gene_id']

    return pd.Series([chrom, start, end, ensID, "."])

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
    parser.add_argument('--protein_coding', help=("Protein coding only"), action='store_true')
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
