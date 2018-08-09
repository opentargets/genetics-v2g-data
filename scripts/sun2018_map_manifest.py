#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Ed Mountjoy (June 2018)

'''

import argparse
import pandas as pd
import sys

def main():

    # Parse args
    args = parse_args()

    #
    # Prepare map data
    #

    # Load map
    map = pd.read_csv(args.map, sep='\t', header=None)
    map.columns = ['gene', 'transcript', 'tss', 'hgnc', 'uniprot', 'chrom']

    # Only keep standard chromosomes
    valid_chroms = set([str(chrom) for chrom in range(1, 23)] +
                        ['X', 'Y', 'MT'])
    map = map.loc[map.chrom.isin(valid_chroms), :]

    # Calculate median TSS for each gene
    gene_tss = (map.loc[:, ['chrom', 'gene', 'tss']]
                   .groupby(['chrom', 'gene'])
                   .median()
                   .reset_index())

    # Assert only 1 gene
    assert(gene_tss.gene.duplicated().sum() == 0)

    # Merge
    # map2 = pd.merge(map.loc[:, ['gene', 'chrom', 'hgnc', 'uniprot']].drop_duplicates(),
    #                 gene_tss,
    #                 how='left')

    #
    # Prepare manifest
    #

    # Load manifest
    mani = pd.read_csv(args.inf, sep=',', header=0)

    # Make HGNC column
    mani['HGNC'] = mani.SOMAMER_ID.apply(lambda x: x.rsplit('.', 3)[0])

    #
    # Map ensembl gene ids to manifest
    #

    # Make map dictionary
    uniprot_map = dict(map.loc[~pd.isnull(map.uniprot), ['uniprot', 'gene']].values.tolist())
    hgnc_map = dict(map.loc[~pd.isnull(map.hgnc), ['hgnc', 'gene']].values.tolist())
    combined_map = {**uniprot_map, **hgnc_map}

    # Map each row
    mani['ensembl_id'] = mani.apply(map_ensembl_to_mani, map=combined_map, axis=1)

    #
    # Prepare output
    #

    # Split rows by ',' in ensembl column
    mani_split = tidy_split(mani, 'ensembl_id', ',')

    # Merge with uniprot id map
    merged = pd.merge(mani_split, gene_tss,
                     left_on='ensembl_id', right_on='gene',
                     how='left')

    # Write
    merged.to_csv(args.outf, sep='\t', index=None)

    return 0

def map_ensembl_to_mani(row, map):
    ''' Maps Ensembl IDs to the manifest file. Maps both using uniprot IDs and
        HGNC IDs.
    Args:
        row (Pandas series): single row from manifest df
        map (dict): Dict mapping both hgnc and uniprot to ensembl
    Returns:
        comma separated list of ensembl IDs
    '''
    ens_ids = []
    # Map on uniprot
    if not pd.isnull(row.UniProt):
        for uniprot in row.UniProt.split(','):
            if uniprot in map:
                ens_ids.append(map[uniprot])
    # Map on hgnc
    if not pd.isnull(row.HGNC):
        for hgnc in row.HGNC.split('.'):
            if hgnc in map:
                ens_ids.append(map[hgnc])
    return ','.join(list(set(ens_ids)))

def tidy_split(df, column, sep='|', keep=False):
    """
    Split the values of a column and expand so the new DataFrame has one split
    value per row. Filters rows where the column is missing.

    Params
    ------
    df : pandas.DataFrame
        dataframe with the column to split and expand
    column : str
        the column to split and expand
    sep : str
        the string used to split the column's values
    keep : bool
        whether to retain the presplit value as it's own row

    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    df = df.dropna(subset=[column])
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('Input file'), type=str, required=True)
    parser.add_argument('--outf', metavar="<file>", help=("Output file"), type=str, required=True)
    parser.add_argument('--map', metavar="<file>", help=("Biomart map file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
