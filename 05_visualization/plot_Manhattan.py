# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 20:28:05 2018

@author: Caiyd
"""

import json
import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import cycle


def load_allchrom_data(infile, chr_col, loc_col, val_col, log_trans, neg_logtrans):
    """
    infile contain several chromosomes
    """
    df = pd.read_csv(infile,
                     sep='\t',
                     usecols = [chr_col, loc_col, val_col],
                     dtype={chr_col: str,
                            loc_col: float,
                            val_col: float})
    if log_trans:
        df.loc[:, val_col] = np.log10(df[val_col].values)
    elif neg_logtrans:
        df.loc[:, val_col] = -np.log10(df[val_col].values)
    df.dropna(inplace=True)
    return df


def plot(df, chr_col, loc_col, val_col, xlabel, ylabel, ylim, invert_yaxis, top_xaxis, cutoff, highlight, outfile, ticklabelsize, figsize, axlabelsize):
    sns.set_style('white', {'ytick.major.size': 3, 'xtick.major.size': 3})
    fig, ax = plt.subplots(1, 1, figsize=figsize)
#    offsets = {}
    loc_offset = 0
    xticks = []
    xticklabels = []
    for chrom, color in zip(df[chr_col].unique(), cycle(['#4C72B0', '#DD8452'])):
#        offsets[chrom] = loc_offset
        tmpdf = df.loc[df[chr_col]==chrom, [loc_col, val_col]]
        tmpdf[loc_col] += loc_offset
        xticklabels.append(chrom)
        xticks.append(tmpdf[loc_col].median())
        tmpdf.plot(kind='scatter', x=loc_col, y=val_col, ax=ax, s=6, color=color, marker='o')
        if isinstance(highlight, pd.DataFrame):
            hdf = highlight.loc[highlight[chr_col]==chrom, [loc_col, val_col]]
            if hdf.shape[0] > 0:
                hdf[loc_col] += loc_offset
                hdf.plot(kind='scatter', x=loc_col, y=val_col, ax=ax, s=6, color='#55A868', marker='o')
        if cutoff:
            ax.hlines(cutoff[chrom], tmpdf[loc_col].values[0], tmpdf[loc_col].values[-1])
        loc_offset = tmpdf[loc_col].values[-1] # assume loc is sorted
    ax.set_xlabel(xlabel, fontsize=axlabelsize)
    ax.set_ylabel(ylabel, fontsize=axlabelsize)
    plt.xticks(xticks, xticklabels)
    ax.set_xlim([df[loc_col].values[0], tmpdf[loc_col].values[-1]])
    if ylim:
        ax.set_ylim(ylim)
    if invert_yaxis:
        ax.invert_yaxis()
    if top_xaxis:
        ax.xaxis.tick_top()
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.xaxis.set_label_position('top')
    else:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(ticklabelsize)
    plt.savefig(f'{outfile}', dpi=300)
    plt.close()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--infile', help='tsv文件,包含header')
@click.option('--chr-col', help='染色体列名')
@click.option('--loc-col', help='x轴值列名')
@click.option('--val-col', help='y轴值列名')
@click.option('--log-trans', is_flag=True, default=False, help='对val列的值取以10为底的对数')
@click.option('--neg-logtrans', is_flag=True, default=False, help='对val列的值取以10为底的负对数(画p值)')
@click.option('--outfile', help='输出文件,根据拓展名判断输出格式')
@click.option('--xlabel', help='输入散点图x轴标签的名称')
@click.option('--ylabel', help='输入散点图y轴标签的名称')
@click.option('--ylim', nargs=2, type=float, default=None, help='y轴的显示范围,如0 1, 默认不限定')
@click.option('--invert-yaxis',  is_flag=True, default=False, help='flag, 翻转y轴, 默认不翻转')
@click.option('--top-xaxis',  is_flag=True, default=False, help='flag, 把x轴置于顶部, 默认在底部')
@click.option('--cutoff', default=None, help='json格式的cutoff, 脚本cal_norm_isf.py的计算结果')
@click.option('--highlight', default=None, help='和infile相同格式的文件,在该文件中的点会在曼哈顿图中单独高亮出来,default=None')
@click.option('--ticklabelsize', help='刻度文字大小', default=10)
@click.option('--figsize', nargs=2, type=float, help='图像长宽, 默认20 5', default=(20, 5))
@click.option('--axlabelsize', help='x轴y轴label文字大小', default=10)
def main(infile, chr_col, loc_col, val_col, log_trans, neg_logtrans, outfile,
         xlabel, ylabel, ylim, invert_yaxis, top_xaxis, cutoff, highlight, ticklabelsize, figsize, axlabelsize):
    """
    \b
    曼哈顿图
    使用infile中的chr_col和loc_col列作为x轴, 对应的val_col值画在y轴
    输出outfile, 根据outfile指定的拓展名输出相应的格式
    """
    print(__doc__)
    print(main.__doc__)
    df = load_allchrom_data(infile, chr_col, loc_col, val_col, log_trans, neg_logtrans)
    if highlight:
        highlight = load_allchrom_data(highlight, chr_col, loc_col, val_col, log_trans, neg_logtrans)
    if cutoff:
        with open(cutoff) as f:
            cutoff = json.load(f)
    plot(df, chr_col, loc_col, val_col, xlabel, ylabel, ylim, invert_yaxis, top_xaxis, cutoff, highlight, outfile, ticklabelsize, figsize, axlabelsize)

if __name__ == '__main__':
    main()
