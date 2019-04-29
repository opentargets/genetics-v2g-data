#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Ed Mountjoy

LiftOver an interval file from one build to another
'''

import sys
import argparse
import gzip
from pyliftover import LiftOver

def main():

    # Parse args
    args = parse_args()
    allowed_chroms = set([str(x) for x in range(1, 23)] + ['X', 'Y', 'MT'])

    # Initate liftover
    lo = LiftOver(args.chainfile)

    # Liftover interval file
    with gzip.open(args.inf, 'r') as in_h:
        with gzip.open(args.outf, 'w') as out_h:

            # Copy header to new file
            out_h.write(in_h.readline())

            # Process each line of the input
            for line in in_h:
                parts = line.decode().rstrip().split('\t')
                
                # LiftOver coordinates
                start_lo = remove_invalid_chromosomes(
                    lo.convert_coordinate(parts[0], int(parts[1]) - 1),
                    allowed_chroms,
                    strip_chr=True)
                end_lo = remove_invalid_chromosomes(
                    lo.convert_coordinate(parts[0], int(parts[2]) - 1),
                    allowed_chroms,
                    strip_chr=True)

                # Ignore line if None returned (input chrom not recognised)
                if start_lo is None or end_lo is None:
                    continue

                # Ignore line if no targets found in liftover
                if len(start_lo) == 0 or len(end_lo) == 0:
                    continue
                
                # Take first target for start_lo
                chrom = start_lo[0][0]
                start_pos = start_lo[0][1] + 1

                # For end_lo take first target that matches `chrom`
                end_pos = None
                for record in end_lo:
                    if record[0] == chrom:
                        end_pos = record[1] + 1
                        break
                
                # Ignore line if no matching end position was found
                if end_pos is None:
                    continue
                
                # Assert that the start < end
                if start_pos >= end_pos:
                    start_pos, end_pos = end_pos, start_pos

                # Check difference in length between original and liftover 
                # intervals
                old_len = abs(int(parts[1]) - int(parts[2]))
                new_len = abs(start_pos - end_pos)
                if abs(old_len - new_len) > args.maxdiff:
                    continue
                
                # Replace co-ordinates
                parts[0] = str(chrom)
                parts[1] = str(start_pos)
                parts[2] = str(end_pos)

                # Write output
                out_h.write(('\t'.join(parts) + '\n').encode())
    
    return 0

def remove_invalid_chromosomes(lo, allowed_chroms, strip_chr=True):
    ''' Process the output of pyliftover to remove invalid chromosomes
    Params:
        lo (list): list of pyliftover results
        allowed_chroms (set): set of allowed chrom names
        strip_chr (bool): whether to remove chr from chromosome names
    Returns:
        list of pyliftover results
    '''
    lo_clean = []
    for record in lo:
        record = list(record)
        if strip_chr:
            record[0] = record[0].replace('chr', '')
        if record[0] in allowed_chroms:
            lo_clean.append(record)
    return lo_clean

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('Input file'), type=str, required=True)
    parser.add_argument('--outf', metavar="<file>", help=("Output file"), type=str, required=True)
    parser.add_argument('--chainfile', metavar="<file>", help=("Chain file for liftover"), type=str, required=True)
    parser.add_argument('--maxdiff', metavar="<int>",
                        help=("Max difference between old and new interval lengths"),
                        type=int, required=False, default=100)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
