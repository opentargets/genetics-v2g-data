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
# from scipy.stats import rankdata

def main():

    # Parse args
    args = parse_args()

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
    data.to_csv('output/test_out_data.tsv', sep='\t', index=None) # DEBUG

    # Print stats
    print('Num unique varids: ', data.variant_id.nunique())

    #
    # Get score per source -----------------------------------------------------
    #

    # Aggregate per source, combining max feature score and non-null feature count
    print('Creating max and count scores per feature...')
    ft_grp_agg = (
                  data.groupby(['variant_id', 'source_id', 'gene_name'])
                      .score_scaled
                      .agg(
                           [np.max,
                            lambda x: np.count_nonzero(~np.isnan(x))]
                          )
                      .rename(columns={'amax':'ft_max', '<lambda>':'ft_count'})
                      .stack()
                      .reset_index()
                      .rename(columns={'level_3':'score', 0:'value'})
                 )

    # Scale rank from 0-1 to allow averaging across score types
    print('Scaling scores 0-1 per score type...')
    ft_grp_agg['value_scale'] = ( ft_grp_agg.groupby(['variant_id', 'source_id', 'score'])
                                           .value
                                           .transform(lambda x: minmax_scale(x.astype(float))) )
    ft_grp_agg.to_csv('output/test_out_agg.tsv', sep='\t') # DEBUG

    # Average over ft_max, ft_count using per source-ft weights
    print('Weighted mean over score types...')
    ft_weights = {
        'vep':           {'ft_max': 1,   'ft_count': 0},
        'javierre2016':  {'ft_max': 0.8, 'ft_count': 0.2},
        'andersson2014': {'ft_max': 0.8, 'ft_count': 0.2},
        'thurman2012':   {'ft_max': 0.8, 'ft_count': 0.2},
        'gtex_v7':       {'ft_max': 0.8, 'ft_count': 0.2}
    }
    ft_grp_agg_overall = (
                ft_grp_agg.groupby(['variant_id', 'source_id', 'gene_name'])
                          .apply(weighted_mean_ft, weight_map=ft_weights)
                          .reset_index()
                          .rename(columns={0:'ft_agg'})
                )
    ft_grp_agg_overall.to_csv('output/test_out_agg_overall.tsv', sep='\t', index=None) # DEBUG

    #
    # Get score per variant_id -------------------------------------------------
    #

    # Average over source_ids using per source weights
    print('Weighted mean over source types...')
    source_weights = {
        'vep': 0.2,
        'javierre2016': 0.2,
        'andersson2014': 0.2,
        'thurman2012': 0.2,
        'gtex_v7': 0.2
    }
    src_grp_agg = (
                ft_grp_agg_overall.groupby(['variant_id', 'gene_name'])
                          .apply(weighted_mean_src, weight_map=source_weights)
                          .reset_index()
                          .rename(columns={0:'src_agg'})
                )
    src_grp_agg.to_csv('output/test_out_src_agg.tsv', sep='\t', index=None) # DEBUG
    print(src_grp_agg.head())



    return 0

def weighted_mean_ft(grp, weight_map):
    ''' Calculate weighted mean using (source_id, score_type) specific weights.
    Args:
        grp (pd df): group from pd.groupby().apply()
        weight_map (dict): dict of weights source_id->score_type->weight
    Returns:
        Weighted mean across score_types
    '''
    weights = [weight_map[source][score] for source, score in zip(grp.source_id, grp.score)]
    weighted_mean = np.average(grp.value_scale, weights=weights)
    return weighted_mean

def weighted_mean_src(grp, weight_map):
    ''' Calculate weighted mean using source_id specific weights.
    Args:
        grp (pd df): group from pd.groupby().apply()
        weight_map (dict): dict of weights source_id->weight
    Returns:
        Weighted mean across score_types
    '''
    weights = [weight_map[source] for source in grp.source_id]
    weighted_mean = np.average(grp.ft_agg, weights=weights)
    return weighted_mean

def quantiles(l):
    ''' Converts a list/series of values to quantiles of same shape.
        sklearn.quantile_transform() requires data to be reshaped
    '''
    return quantile_transform(l.as_matrix().reshape(-1, 1)).reshape(1, -1)[0]

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
    parser.add_argument('--out', metavar="<str>", type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
