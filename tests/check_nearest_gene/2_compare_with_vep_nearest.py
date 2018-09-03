#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import pandas as pd

def main():

    # Args
    in_ot = '../../output/nearest_gene/180830/nearest_gene.tsv.gz'
    in_vep = 'output/vcf_nearest_gene.txt'

    # Load
    ot = pd.read_csv(in_ot, sep='\t', header=0, nrows=1000000)
    vep = pd.read_csv(in_vep, sep='\t', comment='#', header=None).iloc[:, [1, 13]]
    vep.columns = ['loc', 'info']

    # Parse VEP fields
    vep['vep_chrom'], vep['vep_pos'] = vep['loc'].str.split(':').str
    vep['vep_nearest'] = vep['info'].apply(lambda x: x.split('NEAREST=')[-1])
    vep = vep.loc[:, ['vep_chrom', 'vep_pos', 'vep_nearest']]

    # Parse ot fields
    ot['ot_chrom'], ot['ot_pos'], ot['ot_ref'], ot['ot_alt'] = (
        ot.varid.str.split('_').str
    )

    # Merge ot and vep
    merged = pd.merge(ot, vep, left_on=['ot_chrom', 'ot_pos'],
                      right_on=['vep_chrom', 'vep_pos'], how='inner')
    merged = merged.drop_duplicates(subset=['varid'])

    # Find none matching rows
    non_match = (merged.ensemblid_protein_coding != merged.vep_nearest)
    print("There are {0} of {1} variants that don't match the VEP output".format(
        non_match.sum(), merged.shape[0]))

    return 0

if __name__ == '__main__':

    main()
