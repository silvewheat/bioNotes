# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 15:08:58 2016

@author: Caiyd
"""
def load_fa(fafile):
    seqdict = {}
    with open(fafile) as f:
        seq_id = f.readline().strip().split()[0][1:]
        tmp_seq = []
        for line in f:
            if line[0] != '>':
                tmp_seq.append(line.strip())
            else:
                seqdict[seq_id] = ''.join(tmp_seq)
                seq_id = line.strip().split()[0][1:]
                tmp_seq = []
        seqdict[seq_id] = ''.join(tmp_seq)
    return seqdict


def load_contigs(contigsfile):
    return set([x.strip() for x in open(contigsfile).readlines()])

def Out(seq, contigset):
    for k in seq.keys():
        if k in contigset:
            print(f'>{k}\n{seq[k]}\n', end="")

import sys
seq = load_fa(sys.argv[1])
contigset = load_contigs(sys.argv[2])
Out(seq, contigset)
