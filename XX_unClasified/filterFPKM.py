# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 15:05:10 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np
import pandas as pd
import json



def load_group(groupfile):
    with open(groupfile) as f:
        tgroups = json.load(f)
    return tgroups



@click.command()
@click.option('--genefpkm', help='基因表达量表格')
@click.option('--groupjson', help='每个group包含哪些样本')
@click.option('--outprefix', help='输出前缀')
def main(genefpkm, groupjson, outprefix):
    """
    保留在至少一个group中有10%以上的样本FPKM大于等于1
    genefpkm是基因表达量文件
    gene chr start end name sample1 sample2 ...
    """
    tgroups = load_group(groupjson)
    df = pd.read_csv(genefpkm, sep='\t', low_memory=False)
    print(df.shape)
    groups = []
    samples = []
    print('filter...')
    for group, sample in tgroups.items():
        print(group)
        groups.append(group)
        samples.extend(sample)
        df[group] = np.sum(df[sample] >= 1, axis=1) / len(sample) # 在各个group中表达量大于1的百分比
    # 至少一个group中FPKM大于等于1的样本达到10%
    df = df.loc[np.sum(df[groups] >= 0.1, axis=1) != 0, :]
    print('after filter')
    print(df.shape)
    print('write filtered FPKM')
    df[['chr', 'start', 'end', 'name', 'gene', *samples]].to_csv(f'{outprefix}_sampleFPKM.tsv.gz',
                                                                 sep='\t', index=False, compression='gzip')
    # max
    for group, sample in tgroups.items():
        print(f'max {group}')
        df[group] = df[sample].max(axis=1)
    df[['chr', 'start', 'end', 'name', 'gene', *groups]].to_csv(f'{outprefix}_maxFPKM.tsv.gz',
                                                                 sep='\t', index=False, compression='gzip')
    # median
    for group, sample in tgroups.items():
        print(f'median {group}')
        df[group] = df[sample].median(axis=1)
    df[['chr', 'start', 'end', 'name', 'gene', *groups]].to_csv(f'{outprefix}_medianFPKM.tsv.gz',
                                                                 sep='\t', index=False, compression='gzip')
    # mean
    for group, sample in tgroups.items():
        print(f'mean {group}')
        df[group] = df[sample].mean(axis=1)
    df[['chr', 'start', 'end', 'name', 'gene', *groups]].to_csv(f'{outprefix}_meanFPKM.tsv.gz',
                                                                 sep='\t', index=False, compression='gzip')

if __name__ == '__main__':
    main()
