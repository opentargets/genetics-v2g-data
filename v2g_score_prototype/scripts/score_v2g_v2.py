#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import argparse
import pandas as pd
import numpy as np
from sklearn.preprocessing import quantile_transform
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

# from scipy.stats import rankdata

def main():

    # Parse args
    args = parse_args()

    # Source aggregate weights
    source_weights = {
        'vep': 1,
        'javierre2016': 1,
        'andersson2014': 1,
        'thurman2012': 1,
        'gtex_v7': 1,
        'sun2018': 1
    }

    #
    # Load data ----------------------------------------------------------------
    #

    # Load data
    data = load_slice(args.v2g_slice, nrows=None)
    print('Total rows: ', data.shape[0])
    print('size (bp) : ', data.position.max() - data.position.min())

    # Load vep map
    vep_map = load_vep_map()

    #
    # Derive features ----------------------------------------------------------
    #

    print('Preprocessing...')

    # Add score column
    data['score'] = data.apply(create_score_column, vep_map=vep_map, axis=1)

    # Wipe VEP features
    data.loc[data.type_id == 'vep', 'feature'] = 'not_a_feature'

    # Remove rows with null scores
    data = data.loc[~pd.isnull(data.score), :]

    # Calculate quantiles per feature
    data['score_scaled'] = (
        data.groupby(('source_id', 'feature'))
            .score
            .transform(lambda x: quantiles(x, n=10))
        )
    outf = '{0}.preprocessed.tsv.gz'.format(args.outpref)
    data.to_csv(outf, sep='\t', index=None, compression='gzip')

    # Print stats
    print('Num unique varids: ', data.variant_id.nunique())

    #
    # Get score per source -----------------------------------------------------
    #

    # Aggregate per source, calculate max feature score
    print('Creating max and count scores per feature...')
    per_source_scores = (
        data.groupby(['variant_id', 'source_id', 'gene_name'])
            .score_scaled
            .agg([np.max])
            .rename(columns={'amax':'ft_max'})
        )

    # Melt into long format
    per_source_scores = (
        per_source_scores.stack()
                         .reset_index()
                         .rename(columns={'level_3':'score_type', 0:'score'})
        )
    outf = '{0}.per_source_scores.tsv'.format(args.outpref)
    per_source_scores.to_csv(outf, sep='\t', index=None)

    #
    # Get score per variant_id -------------------------------------------------
    #

    # Average over source_ids using per source weights
    print('Weighted mean over source types...')
    per_varid = (
        per_source_scores.groupby(['variant_id', 'gene_name'])
                         .apply(weighted_mean_src, score_col='score', weight_map=source_weights)
                         .reset_index()
                         .rename(columns={0:'src_agg'})
        )
    outf = '{0}.per_variant_scores.tsv'.format(args.outpref)
    per_varid.to_csv(outf, sep='\t', index=None)

    #
    # With variant fixed, make heatmaps ----------------------------------------
    #

    # Set matplotlib to auto layout
    rcParams.update({'figure.autolayout': True})

    # args.varid = '16_53563398_C_A'

    # Stop here if varid is not specified
    if not args.varid:
        sys.exit()
    elif not args.varid in data.variant_id.unique():
        sys.exit('--varid {0} not found in data'.format(args.varid))

    # Make one heatmap per source_id
    for name, grp in data.loc[data.variant_id == args.varid, :].groupby('source_id'):
        plot_data = grp.pivot_table(index='gene_name',
                                    columns='feature',
                                    values='score_scaled')
        ax = sns.heatmap(plot_data, cmap=sns.color_palette("Blues"), annot=True)
        ax.set_title('{0}\n{1}'.format(args.varid, name))
        out_plot = '{0}.plot.{1}.per_tissue.png'.format(args.outpref, name)
        ax.figure.savefig(out_plot)
        plt.clf()

    # Make heatmap summarising source_ids, using aggregate score
    plot_data = (
        per_source_scores.loc[(per_source_scores.variant_id == args.varid), :]
                         .pivot_table(index='gene_name',
                                      columns='source_id',)
                                      # values='score')
        )
    ax = sns.heatmap(plot_data, cmap=sns.color_palette("Blues"), annot=True)
    ax.set_title('{0}\n{1}\n{2}'.format(args.varid, 'source_id_summary', 'aggregate_score'))
    out_plot = '{0}.plot.per_source.png'.format(args.outpref)
    ax.figure.savefig(out_plot)
    plt.clf()

    return 0

def weighted_mean_src(grp, score_col, weight_map):
    ''' Calculate weighted mean using source_id specific weights.
    Args:
        grp (pd df): group from pd.groupby().apply()
        weight_map (dict): dict of weights source_id->weight
    Returns:
        Weighted mean across score_types
    '''
    weights = [weight_map[source] for source in grp.source_id]
    weighted_mean = np.average(grp[score_col], weights=weights)
    return weighted_mean

def quantiles(l, n=10):
    ''' Converts a list/series of values to quantiles of same shape.
        sklearn.quantile_transform() requires data to be reshaped
    '''
    qt = quantile_transform(X=l.as_matrix().reshape(-1, 1),
                            n_quantiles=n)
    return qt.reshape(1, -1)[0]

def create_score_column(row, vep_map):
    ''' Returns a score specific to type_id
    '''
    if row.type_id == 'pchic':
        if row.interval_score != 0:
            return row.interval_score
        else:
            return np.nan
    elif row.type_id == 'fantom5':
        if row.interval_score != 0:
            return row.interval_score
        else:
            return np.nan
    elif row.type_id == 'dhscor':
        if row.interval_score != 0:
            return row.interval_score
        else:
            return np.nan
    elif row.type_id == 'vep':
        return vep_map[row.feature]
    elif row.type_id == 'eqtl':
        return 1 - row.qtl_pval
    else:
        sys.exit('Error: type_id not recognised {0}'.format(row.type_id))

def load_slice(inf, nrows=None):
    '''
    Args:
        inf (str): file containing slice of v2g
    Returns:
        pd.df
    '''
    df = pd.read_csv(inf, sep='\t', header=None, nrows=nrows)
    df.columns = [
        'chr_id',
        'position',
        'ref_allele',
        'alt_allele',
        'variant_id',
        'rs_id',
        'gene_chr',
        'gene_id',
        'gene_start',
        'gene_end',
        'gene_name',
        'feature',
        'type_id',
        'source_id',
        'csq_counts',
        'qtl_beta',
        'qtl_se',
        'qtl_pval',
        'interval_score']
    return df

def load_vep_map():
    vm_df = pd.read_csv('docs/vep_consequences_180807.tsv', sep='\t', header=0)
    vep_map = dict(zip(vm_df['SO term'], vm_df['v2g_score']))
    return vep_map

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--v2g_slice', metavar="<file>", type=str, required=True)
    parser.add_argument('--varid', metavar="<file>", type=str)
    parser.add_argument('--outpref', metavar="<str>", help='Output prefix', type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
