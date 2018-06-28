# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 15:56:57 2018

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
@click.option('--grouplist', help='所有group的名字')
@click.option('--targetgroup', help='找这个group富集的基因')
@click.option('--fold', help='多少倍算富集', type=int)
@click.option('--maxdiff', help='例外group上限', type=int)
@click.option('--outprefix', help='输出前缀')
def main(genefpkm, grouplist, targetgroup, fold, maxdiff, outprefix):
    """
    保留在至少一个group中有10%以上的样本FPKM大于等于1
    genefpkm是基因表达量文件
    gene chr start end name sample1 sample2 ...
    """
    df = pd.read_csv(genefpkm, sep='\t', low_memory=False)
    groups = [x.strip() for x in open(grouplist).readlines()]
    othergroups = groups.copy()
    othergroups.remove(targetgroup)
    target_exp_fold = df[targetgroup] / fold
    # 不满足条件的group数量
    x = target_exp_fold.shape[0]
    diffcount = np.sum(df[othergroups].values >= target_exp_fold.values.reshape(x, 1), axis=1)
    df.loc[(diffcount <= maxdiff) & (df[targetgroup] >= 1), :].to_csv(f'{outprefix}.tsv.gz', sep='\t', index=False, compression='gzip')


if __name__ == '__main__':
    main()

