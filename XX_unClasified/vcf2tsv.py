# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 10:19:44 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import os
import click
import numpy as np
import pandas as pd


def load_bcf2df(bcffile, region=''):
    snp_list = []
    index_list = []
    cmd = 'bcftools view %s %s' % (bcffile, region)
    for line in os.popen(cmd):
        if line[0] != '#':
            line = line.strip().split()
            snp_list.append([x.split(':')[0] for x in line[9:]])
            index_list.append('%s:%s' % (line[0], line[1]))
        elif line[:6] == '#CHROM':
            name_list = line.strip().split('\t')[9:]
    return pd.DataFrame(snp_list, columns=name_list, index=index_list)


def load_sample_order(order_file):
    sample_order_list = []
    with open(order_file) as f:
        for line in f:
            if line[0] != '#':
                line = line.strip()
                sample_order_list.append(line)
    return sample_order_list




@click.command()
@click.option('--bcffile', help='输入的bcf文件')
@click.option('--region', help='要提取的区域，如12:1000-2000')
@click.option('--orderfile', help='个体ID顺序，输出文件会按这个排')
@click.option('--outfile', help='输出文件名')
def main(bcffile, region, orderfile, outfile):
    """
    把vcf文件中指定区域转成tsv
    """
    df = load_bcf2df(bcffile, region)
    sample_order_list = load_sample_order(orderfile)
    df[sample_order_list].to_csv(outfile, sep='\t')



if __name__ == '__main__':
    main()
