#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Ed Mountjoy (June 2018)

Converts Javierre et al file to bed format
'''

import argparse
import pandas as pd

def main():

    # Parse args
    args = parse_args()

    # Load
    df = pd.read_csv(args.inf, sep='\t', header=None)

    # Split column 3 into separate cols
    expanded = df.iloc[:, 3].str.split(':|-|,', expand=True)
    df = pd.concat([df.iloc[:, 0:3], expanded, df.iloc[:, 4:6]], axis=1)
    df.columns = range(df.shape[1])

    # Remove 'chr' from 0 and 3
    for col in [0, 3]:
        df.iloc[:, col] = df.iloc[:, col].str.replace('chr', '')

    # Save
    df.to_csv(args.outf, sep='\t', index=None, header=None, compression='gzip')

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
