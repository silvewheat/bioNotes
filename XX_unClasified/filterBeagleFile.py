# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 11:15:12 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import sys
import gzip
import click
import numpy as np
import pandas as pd


def print_chunk(n):
    sys.stdout.write(' ' * 30 + '\r')
    sys.stdout.flush()
    sys.stdout.write(f'chunk {n}' + '\r')
    sys.stdout.flush()


@click.command()
@click.option('--beaglefile', help='压缩后的beagle文件')
@click.option('--sites', help='只保留出现在该文件中的位点,与beagle中的marker一致,(可选)', default=None)
@click.option('--samples', help='只保留出现在该文件中的个体,(可选)', default=None)
@click.option('--rmtansi', help='去掉所有转换(transition), (flag)', is_flag=True, default=False)
@click.option('--chuncksize', help='一次处理XX行，默认10000', default=10000, type=int)
@click.option('--outfile', help='输出压缩后文件的名字(xx.gz)')
def main(beaglefile, sites, samples, rmtansi, chuncksize, outfile):
    """
    从bealge.gz文件中提取指定个体指定位点
    beagle格式：http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods#Beagle_format
    allele codes as 0=A, 1=C, 2=G, 3=T
    """
    if samples:
        querysamples = [x.strip() for x in open(samples).readlines()]
    if sites:
        sites = set(pd.read_csv(sites, sep='\t', usecols=[0], header=None).T.values[0].tolist())
    usecols = [0, 1, 2]
    with gzip.open(beaglefile, 'rb') as f:
        header = np.array(f.readline().decode().split())
        index = np.array(list(range(len(header))))
        if samples:
            for smid in querysamples:
                usecols.extend(index[header == smid])
        else:
            usecols = list(index)
    with gzip.open(outfile, 'wb') as f:
        oheader = '\t'.join(header[usecols]) + '\n'
        f.write(oheader.encode())
    reader = pd.read_csv(beaglefile, sep='\t', skiprows=1, header=None, iterator=True)
    loop = True
    nchunk = 0
    while loop:
        nchunk += 1
        print_chunk(nchunk)
        try:
            chunk = reader.get_chunk(chuncksize)
            if sites:
                chunk = chunk.loc[chunk[0].isin(sites), :]
            if samples:
                chunk = chunk.loc[:, usecols]
            if rmtansi:
                chunk = chunk.loc[np.fabs(chunk[1] - chunk[2])==2, :]
            chunk.to_csv(outfile, sep='\t', index=False, header=False, compression='gzip', mode='ab')
        except StopIteration:
            loop = False
            print('Done!')


if __name__ == '__main__':
    main()

