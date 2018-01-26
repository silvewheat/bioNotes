# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 11:18:47 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np
import pandas as pd
from numpy.linalg import eig



@click.command()
@click.option('--covfile', help='pcAngsd输出的.cov文件')
@click.option('--k', help='输出前k个主成分, 默认前5个', default=5)
@click.option('--outprefix', help='输出文件前缀')
def main(covfile, k, outprefix):
    """
    eigendecomposition
    """
    cov = np.genfromtxt(covfile)
    eigenvalues, eigenvectors = eig(cov)
    sumval = np.sum(eigenvalues)
    eigenvec = eigenvectors.T[eigenvalues.argsort()[::-1]] # 要转秩 每列是一个vec
    eigenval = np.sort(eigenvalues)[::-1][:k]
    pcs = {}
    for i in range(k):
        exp = eigenval[i] / sumval
        pcs[f'PC{i+1}({exp*100:.2f}%)'] = eigenvec[i]
    pd.DataFrame(pcs).to_csv(f'{outprefix}.tsv', index=False, sep='\t')


if __name__ == '__main__':
    main()
