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




def Out(seq):
    for k in sorted(list(seq.keys())):
        length = len(seq[k])
        print(k+'\t'+str(length))

import sys
seq = load_fa(sys.argv[1])
Out(seq)
