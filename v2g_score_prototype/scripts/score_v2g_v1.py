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
from sklearn.preprocessing import minmax_scale
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

# from scipy.stats import rankdata

def main():

    # Parse args
    args = parse_args()

    # Tissue (feature) aggregate max/count weights
    ft_weights = {
        'vep':           {'ft_max': 1,   'ft_count_scaled': 0},
        'javierre2016':  {'ft_max': 0.8, 'ft_count_scaled': 0.2},
        'andersson2014': {'ft_max': 0.8, 'ft_count_scaled': 0.2},
        'thurman2012':   {'ft_max': 0.8, 'ft_count_scaled': 0.2},
        'gtex_v7':       {'ft_max': 0.8, 'ft_count_scaled': 0.2},
        'sun2018':       {'ft_max': 0.8, 'ft_count_scaled': 0.2}
    }

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
            .transform(lambda x: quantiles(x))
        )
    outf = '{0}.preprocessed.tsv.gz'.format(args.outpref)
    data.to_csv(outf, sep='\t', index=None, compression='gzip')

    # Print stats
    print('Num unique varids: ', data.variant_id.nunique())

    #
    # Get score per source -----------------------------------------------------
    #

    # Aggregate per source, calculate max feature score and non-null feature count
    print('Creating max and count scores per feature...')
    per_source_scores = (
        data.groupby(['variant_id', 'source_id', 'gene_name'])
            .score_scaled
            .agg(
                 [np.max,
                  lambda x: np.count_nonzero(~np.isnan(x))]
                )
            .rename(columns={'amax':'ft_max', '<lambda>':'ft_count'})
        )


    # Scale tissue counts from 0-1 to allow averaging with quantiles
    per_source_scores['ft_count_scaled'] = (
        per_source_scores.groupby(['variant_id', 'source_id'])
                         .ft_count
                         .transform(lambda x: minmax_scale(x.astype(float)))
                          )
    # Drop unscaled counts column and melt into long format
    per_source_scores = (
        per_source_scores.drop(labels=['ft_count'], axis=1)
                         .stack()
                         .reset_index()
                         .rename(columns={'level_3':'score_type', 0:'value'})
        )
    outf = '{0}.per_source_scores.tsv'.format(args.outpref)
    per_source_scores.to_csv(outf, sep='\t', index=None)

    # Average over ft_max, ft_count using per source-ft weights
    print('Weighted mean over score types...')
    per_source_agg = (
        per_source_scores.groupby(['variant_id', 'source_id', 'gene_name'])
                         .apply(weighted_mean_ft, score_col='value', weight_map=ft_weights)
                         .reset_index()
                         .rename(columns={0:'ft_agg'})
        )
    outf = '{0}.per_source_aggregate_scores.tsv'.format(args.outpref)
    per_source_agg.to_csv(outf, sep='\t', index=None)

    #
    # Get score per variant_id -------------------------------------------------
    #

    # Average over source_ids using per source weights
    print('Weighted mean over source types...')
    per_varid = (
        per_source_agg.groupby(['variant_id', 'gene_name'])
                      .apply(weighted_mean_src, score_col='ft_agg', weight_map=source_weights)
                      .reset_index()
                      .rename(columns={0:'src_agg'})
        )
    outf = '{0}.per_variant_aggregate_scores.tsv'.format(args.outpref)
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

    # Make heatmap summarising source_ids, using feature max and count
    for score_type in ['ft_max', 'ft_count_scaled']:
        plot_data = (
            per_source_scores.loc[(per_source_scores.variant_id == args.varid) & (per_source_scores.score_type == score_type), :]
                             .pivot_table(index='gene_name',
                                          columns='source_id',
                                          values='value')
            )
        ax = sns.heatmap(plot_data, cmap=sns.color_palette("Blues"), annot=True)
        print(plot_data)
        ax.set_title('{0}\n{1}\n{2}'.format(args.varid, 'source_id_summary', score_type))
        out_plot = '{0}.plot.per_source.{1}.png'.format(args.outpref, score_type)
        ax.figure.savefig(out_plot)
        plt.clf()

    # Make heatmap summarising source_ids, using aggregate score
    plot_data = (
        per_source_agg.loc[(per_source_agg.variant_id == args.varid), :]
                         .pivot_table(index='gene_name',
                                      columns='source_id',
                                      values='ft_agg')
        )
    ax = sns.heatmap(plot_data, cmap=sns.color_palette("Blues"), annot=True)
    ax.set_title('{0}\n{1}\n{2}'.format(args.varid, 'source_id_summary', 'aggregate_score'))
    out_plot = '{0}.plot.per_source.aggregate.png'.format(args.outpref, score_type)
    ax.figure.savefig(out_plot)
    plt.clf()

    return 0

def weighted_mean_ft(grp, score_col, weight_map):
    ''' Calculate weighted mean using (source_id, score_type) specific weights.
    Args:
        grp (pd df): group from pd.groupby().apply()
        weight_map (dict): dict of weights source_id->score_type->weight
    Returns:
        Weighted mean across score_types
    '''
    weights = [weight_map[source][score] for source, score in zip(grp.source_id, grp.score_type)]
    weighted_mean = np.average(grp[score_col], weights=weights)
    return weighted_mean

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

def quantiles(l):
    ''' Converts a list/series of values to quantiles of same shape.
        sklearn.quantile_transform() requires data to be reshaped
    '''
    qt = quantile_transform(X=l.as_matrix().reshape(-1, 1),
                            n_quantiles=100)
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
    vep_map = {
        "cttv_mapping_pipeline": 1.,
        "curator_inference": 1.,
        "trinucleotide_repeat_expansion": 1.,
        "sequence_variant": 0.5,
        "regulatory_region_variant": 0.6,
        "stop_retained_variant": 0.65,
        "splice_acceptor_variant": 0.95,
        "splice_donor_variant": 0.95,
        "stop_lost": 0.9,
        "coding_sequence_variant": 0.95,
        "initiator_codon_variant": 0.7,
        "missense_variant": 0.7,
        "stop_gained": 0.95,
        "frameshift_variant": 0.95,
        "non_coding_transcript_variant": 0.65,
        "mature_miRNA_variant": 0.65,
        "NMD_transcript_variant": 0.65,
        "5_prime_UTR_variant": 0.65,
        "3_prime_UTR_variant": 0.65,
        "incomplete_terminal_codon_variant": 0.9,
        "intron_variant": 0.65,
        "intergenic_variant": 0.5,
        "splice_region_variant": 0.95,
        "upstream_gene_variant": 0.6,
        "downstream_gene_variant": 0.6,
        "TF_binding_site_variant": 0.6,
        "non_coding_transcript_exon_variant": 0.65,
        "protein_altering_variant": 0.7,
        "synonymous_variant": 0.65,
        "inframe_insertion": 0.7,
        "inframe_deletion": 0.7,
        "conservative_inframe_deletion": 0.5,
        "transcript_amplification": 0.6,
        "regulatory_region_amplification": 0.6,
        "TFBS_amplification": 0.6,
        "transcript_ablation": 1.,
        "regulatory_region_ablation": 0.6,
        "TFBS_ablation": 0.6,
        "feature_truncation": 0.6,
        "feature_elongation": 0.6,
        "start_lost": 0.95,
        "Nearest gene counting from 5' end": 0.5,
        "Regulatory nearest gene 5' end": 0.6
    }
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
