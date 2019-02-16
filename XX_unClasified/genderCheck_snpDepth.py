# -*- coding: utf-8 -*-
"""
Created on Sun Sep 30 17:27:35 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""


import pysam
import click
import numpy as np
import pandas as pd
from pandas.compat import StringIO
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture



def samdepth(samfile, sitesfile):
    depth = pysam.depth(samfile, '-aa', '-b', sitesfile)
    depth = pd.read_csv(StringIO(depth), sep='\t', header=None, usecols=[2], squeeze=True)
    return depth.values






@click.group()
def main():
    """
    gender distinguish according to read depth in snp
    """



@main.command('indivstat', short_help='stat indiv depth')
@click.option('--samfile', help='input sam/bam file')
@click.option('--autosites', help='autosome snp site file, chrom\\tloc, suport gziped')
@click.option('--chrxsites', help='chrX snp site file, chrom\\tloc, suport gziped')
@click.option('--chrysites', help='chrY snp site file, chrom\\tloc, suport gziped')
@click.option('--outprefix', help='prefix of output file')
@click.option('--savedp', is_flag=True, default=False, help='flag, save the depth per snp site, default is False')
def indivstat(samfile, autosites, chrxsites, chrysites, outprefix, savedp):
    alnfile = pysam.AlignmentFile(samfile)
    smid = alnfile.header['RG'][0]['SM']
    autodp = samdepth(samfile, autosites)
    print('auto sites: {autodp.shape[0]}')
    # 用众数作为标准化的分母
    mode, count = stats.mode(autodp)
    mode = mode[0]
    count = count[0]
    print(f'auto mode: {mode}. count: {count}')
    if mode > 1:
        xdp = samdepth(samfile, chrxsites)
        print(f'x sites: {xdp.shape[0]}')
        xmedian = np.median(xdp)
        print(f'X median: {xmedian}')

        ydp = samdepth(samfile, chrysites)
        print(f'y sites: {ydp.shape[0]}')
        ymedian = np.median(ydp)
        print(f'Y median: {ymedian}')

        with open(f'{outprefix}.out', 'w') as f:
            f.write(f'{smid}\t{xmedian/mode}\t{ymedian/mode}\n')
    else: # mode == 0, sequence depth is too low, use mean
        automean = np.mean(autodp)
        print('auto depth low than 1, use mean')
        print(f'auto mean: {automean}.')
        xdp = samdepth(samfile, chrxsites)
        print(f'x sites: {xdp.shape[0]}')
        xmean = np.mean(xdp)
        print(f'X mean: {xmean}')
        ydp = samdepth(samfile, chrysites)
        print(f'y sites: {ydp.shape[0]}')
        ymean = np.mean(ydp)
        print(f'Y mean: {ymean}')

        with open(f'{outprefix}.out', 'w') as f:
            f.write(f'{smid}\t{xmean/automean}\t{ymean/automean}\n')
    if savedp:
        np.save(f'{outprefix}_auto', autodp)
        np.save(f'{outprefix}_x', xdp)
        np.save(f'{outprefix}_y', ydp)



@main.command('popcluster', short_help='stat indiv depth')
@click.option('--statfile', help='merged all indivstat out file. (cat *.out > merge.out)')
@click.option('--outprefix', help='output file prefix')
def popcluster(statfile, outprefix):
    df = pd.read_csv(statfile, sep='\t', header=None, names=['sample', 'X-rate', 'Y-rate'])
    X = df[['X-rate', 'Y-rate']]
    gmm = GaussianMixture(n_components=2).fit(X)
    labels = gmm.predict(X)
    df['label'] = labels
    y_rates = []
    for label in range(2):
        y_rates.append(df.loc[df['label']==label, 'Y-rate'].median())
    male = np.argmax(y_rates)
    df['gender'] = 'female'
    df.loc[df['label']==male, 'gender'] = 'male'
    df['color'] = df['gender'].map({'male': 'blue', 'female': 'red'})
    df.to_csv(f'{outprefix}.tsv', sep='\t', index=False)

    fig, ax = plt.subplots(1,1, figsize=(5,5))
    ax.scatter(df['X-rate'], df['Y-rate'], color=df['color'])
    ax.set_xlabel('X-rate')
    ax.set_ylabel('Y-rate')
    plt.savefig(f'{outprefix}.jpg', dpi=300)






if __name__ == '__main__':
    main()
