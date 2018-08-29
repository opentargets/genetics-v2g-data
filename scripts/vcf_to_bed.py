#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Ed Mountjoy (August 2018)
'''

import argparse
import sys
import gzip

def main():

    # Parse args
    args = parse_args()

    with gzip.open(args.invcf, 'r') as in_vcf, gzip.open(args.outbed, 'w') as out_bed:
        for line in in_vcf:
            line = line.decode().rstrip()
            # Skip headers
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            chrom, pos, _, ref, multi_alt, *_ = parts
            # Spilt multiallelic
            for alt in multi_alt.split(','):
                varid = '_'.join([chrom, pos, ref, alt])
                end = int(pos)
                start = end - 1
                out_row = [chrom, str(start), str(end), varid, '.']
                out_bed.write(('\t'.join(out_row) + '\n').encode())

    print('COMPLETE')

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--invcf', metavar="<file>", help=('Input file'), type=str, required=True)
    parser.add_argument('--outbed', metavar="<file>", help=("Output file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
