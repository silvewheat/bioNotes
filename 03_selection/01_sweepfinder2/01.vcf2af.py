# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 12:42:26 2017

@author: Caiyd
"""

import click
import numpy as np
from pysam import VariantFile



@click.command()
@click.option('--varfile', help='vcf(.gz)或bcf文件')
@click.option('--keep', help='只计算该文件中列出来的个体')
@click.option('--outprefix', help='输出allele frequency文件的前缀')
def main(varfile, keep, outprefix):
    """
    输出sweepfinder2的allele frequency file文件
    要求输入的vcf文件没有缺失
    输出文件中会把alt allele count为0的过滤掉
    """
    varin = VariantFile(varfile)
    samples = [x.strip() for x in open(keep).readlines()]
    varin.subset_samples(samples)
    print(f'keep samples:\n{samples}')
    ss = len(samples) * 2 # sample size
    with open(f'{outprefix}.SF', 'w') as f:
        f.write('position\tx\tn\tfolded\n')
        for rec in varin.fetch():
            gts = [s['GT'] for s in rec.samples.values()]
            gts = np.array(gts, dtype='int8').flatten()
            ac = np.sum(gts) # alt allele count
            if ac > 0:
                f.write(f'{rec.pos}\t{ac}\t{ss}\t1\n')


if __name__ == '__main__':
    main()

