# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 20:33:14 2018

@author: Caiyd
"""

import click
import numpy as np
from scipy.stats import norm
from collections import defaultdict


def load_parameter(stafile):
    paras = {}
    with open(stafile) as f:
        for nline, line in enumerate(f, 1):
            if nline >= 8:
                tline = line.strip().split()
                paras[tline[1]] = float(tline[2])
    return paras


def get_cutoff(rvs: list):
    """
    输入的rvs是符合正态分布的一组随机变量(random variables)
    Percent point function (inverse of cdf — percentiles).
    """
    mean = np.mean(rvs)
    std = np.std(rvs)
    lower_bound = norm.isf(0.005, loc=mean, scale=std)
    upper_bound = norm.isf(0.975, loc=mean, scale=std)
    mid = norm.isf(0.5, loc=mean, scale=std)
    return lower_bound, upper_bound, mid


@click.command()
@click.option('--stafilelist', help='03.collect.pl.sta文件列表')
@click.option('--outfile', help='输出结果')
def main(stafilelist, outfile):
    all_paras = defaultdict(list)
    with open(stafilelist) as f:
        for file in f:
            paras = load_parameter(file.strip())
            for k, v in paras.items():
                all_paras[k].append(v)
    with open(outfile, 'w') as f:
        f.write('item\t0.025\t0.975\t0.5\n')
        for k, v in all_paras.items():
            lower_bound, upper_bound, mid = get_cutoff(v)
            f.write(f'{k}\t{lower_bound}\t{upper_bound}\t{mid}\n')


if __name__ == '__main__':
    main()





