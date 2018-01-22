# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 21:39:21 2018

@author: Caiyd
"""

import re
import click
import numpy as np
import pandas as pd


def load_annovar(infile):
    """
    annovar, variant_function file
    """
    df = pd.read_csv(infile, sep='\t', header=None,
                     names=['vartype', 'genes', 'chrom', 'start', 'end'],
                     usecols=[0, 1, 2, 3, 4],
                     dtype={'vartype': str,
                            'genes': str,
                            'chrom': str,
                            'start': int,
                            'end': int})
    return df


@click.command()
@click.option('--annfile', help='annovar的注释结果')
@click.option('--queryfile', help='bed格式的注释结果')
@click.option('--noheader', is_flag=True, default=False, help='flag, 默认有header, 没的话加上这个flag')
@click.option('--outfile', help='输出文件名')
def main(annfile, queryfile, noheader, outfile):
    """
    从annovar结果中输出query区域中的基因(chrom, start, end三列一致)
    """
    replace = re.compile(r'\(.*?\)') # 去掉形如FOXP2(XM_018046611.1:c.-279574_-279430delins0),FOXP2(XM_018046611.1:c.-279574_-279430delins0)中的括号
#    excludetypes = {'intergenic', 'upstream', 'downstream', 'upstream;downstream'}
    excludetypes = {'intergenic'}
    ann = load_annovar(annfile)
    ann = ann.loc[~ann['vartype'].isin(excludetypes), :]
    if noheader:
        query = pd.read_csv(queryfile, sep='\t', header=None, usecols=[0,1,2], names=['chrom', 'start', 'end'], dtype={'chrom': str, 'start': int, 'end': int})
    else:
        query = pd.read_csv(queryfile, sep='\t', header=0, usecols=[0,1,2], names=['chrom', 'start', 'end'], dtype={'chrom': str, 'start': int, 'end': int})
    raw_genes = pd.merge(ann, query, how='inner', on=['chrom', 'start', 'end'])['genes']
    genes = []
    for gene in raw_genes:
        genes.extend(replace.sub('', gene).split(','))
    with open(outfile, 'w') as f:
        f.write('\n'.join(np.unique(genes)))

if __name__ == '__main__':
    main()
