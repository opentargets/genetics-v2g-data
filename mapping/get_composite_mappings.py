#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Jeremy Schwartzentruber
#
import pandas as pd
import argparse

def main():
    args = parse_args()
    df = pd.read_json(args.inf, lines=True)
    # Filter df to ensure study is not null
    df = df[df['study'].notnull()]
    df['original_source'] = 'composite_source_feature_hack'
    df['biofeature_string'] = df['study'] + '-' + df['biofeature_string']
    df['biofeature_code'] = df['study'] + '-' + df['biofeature_code']
    df.to_json(args.out, orient='records', lines=True)


def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('Input biofeature mappings'), type=str, required=True)
    parser.add_argument('--out', metavar="<file>", help=('Output biofeature mappings'), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()

df = pd.read_json('eqtl_catalogue.mappings.json', lines=True)