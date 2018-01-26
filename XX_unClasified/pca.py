# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 13:25:00 2017

@author: Caiyd
"""

import click
import pandas as pd
import numpy as np
from numpy.linalg import eig
#from scipy.stats.mstats import zscore



def load_data(infile):
    """
    读纯数字的输入矩阵
    每行一个个体，每列一个特征
    """
    df = pd.read_csv(infile, sep='\t', header=None)
    return df.values


def pca(X, k, standardize):
    if standardize:
        X = X - X.mean(axis = 0)
    X_cov = np.cov(X.T, ddof = 0)
    eigenvalues, eigenvectors = eig(X_cov)
    klarge_index = eigenvalues.argsort()[-k:][::-1] # top k
    k_eigenvectors = eigenvectors[klarge_index]
    return np.dot(X, k_eigenvectors.T), np.sort(eigenvalues)[::-1]


@click.command()
@click.argument('infile')
@click.argument('outprefix')
@click.option('-k', type=int, default=5, help='输出前k个主成分, 默认5个')
@click.option('--standardize/--no-standardize', default=True, help='是否做标准化，默认做')
def main(infile, outprefix, k, standardize):
    """
    ERROR！！！！！！
    输入的INFILE是纯数字矩阵，tab分割，每行一个个体，每列一个特征
    """
    X = load_data(infile)
    pcs, eigenvalues = pca(X, k, standardize)
    pd.DataFrame(pcs).to_csv(f'{outprefix}.pcs', sep='\t', index=False, header=False)
    eigenvalues.tofile(f'{outprefix}.eigval', sep='\n')


if __name__ == '__main__':
    main()
