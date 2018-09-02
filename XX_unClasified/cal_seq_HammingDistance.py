# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 16:21:42 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""


import click
import numpy as np
import pandas as pd
from collections import defaultdict
from itertools import combinations



def cal_dist_excludeN(seq1: np.ndarray, seq2: np.ndarray):
    mask = (seq1 != 'N') & (seq2 != 'N')
    return np.sum(seq1[mask] != seq2[mask])

def cal_dist_all(seq1: np.ndarray, seq2: np.ndarray):
    return np.sum(seq1 != seq2)



def loadfa(fafile):
    seqdict = defaultdict(list)
    with open(fafile) as f:
        seq_id = f.readline().strip().split()[0][1:]
        tmp_seq = []
        for line in f:
            if line[0] != '>':
                tmp_seq.append(list(line.strip()))
            else:
                seqdict[seq_id] = np.concatenate(tmp_seq)
                seq_id = line.strip().split()[0][1:]
                tmp_seq = []
        seqdict[seq_id] = np.concatenate(tmp_seq)
    return seqdict



@click.command()
@click.option('--fafile', help='fasta file, seqence must be in same length!')
@click.option('--rmgap',  is_flag=True, default=False, help='flag, 加上在算距离前会去掉N')
@click.option('--outfile', help='输出文件名')
def main(fafile, rmgap, outfile):
    """
    \b
    cal hamming distance
    seqence must be in same length
    """
    seqdict = loadfa(fafile)
    if rmgap:
        cal_dist = cal_dist_excludeN
    else:
        cal_dist = cal_dist_all
    result = []
    for seqid1, seqid2 in combinations(seqdict.keys(), 2):
        seq1 = seqdict[seqid1]
        seq2 = seqdict[seqid2]
        dist = cal_dist(seq1, seq2)
        result.append([seqid1, seqid2, dist])
        result.append([seqid2, seqid1, dist])
    # diagonal
    for seqid in seqdict.keys():
        result.append([seqid, seqid, 0])
    result = pd.DataFrame(result, columns=['seq1', 'seq2', 'dist'])
    result.pivot_table(values='dist', index='seq1', columns='seq2').to_csv(outfile, sep='\t')


if __name__ == '__main__':
    main()






