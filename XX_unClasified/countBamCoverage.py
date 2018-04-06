# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 14:53:10 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import sys
import click
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
import pysam


def load_region(regionfile):
    regions = []
    with open(regionfile) as f:
        for line in f:
            if line[0] != '#':
                sline = line.strip().split()
                regions.append((sline[0], int(sline[1]), int(sline[2])))
    return regions


def load_alnfile(alnfile, reffile):
    if alnfile.split('.')[-1] == 'cram':
        if not reffile:
            sys.exit("cram file need refrence fasta in --reffile")
        aln = pysam.AlignmentFile(alnfile, reference_filename=reffile)
    else:
        aln = pysam.AlignmentFile(alnfile)
    smid = aln.header['RG'][0]['SM']
    return smid, aln


def count_cov(aln, contig: str, start: int, stop: int):
    """
    start (0-based inclusive)
    stop (0-based inclusive)
    """
    return np.sum(aln.count_coverage(contig=contig, start=start, stop=stop), axis=0)



@click.command()
@click.option('--alnfile', help='cram/bam/sam的文件')
@click.option('--reffile', help='输入的参考基因文件,cram文件必须要,bam不用', default=None)
@click.option('--regionfile', help='需要计算的区域，bed格式')
@click.option('--cutoff', '-t', help='depth coverage, multiple, 默认1, 3, 5, 7, 10', multiple=True, default=None, type=int)
@click.option('--outprefix', help='输出结果文件的前缀')
def main(alnfile, reffile, regionfile, cutoff, outprefix):
    """
    \b
    计算regionfile中的coverage

    @Email: yudongcai216@gmail.com
    """
    cutoffs = (1, 3, 5, 7, 10) if not cutoff else cutoff
    print(f'cutoff: {cutoffs}')
    regions = load_region(regionfile)
    smid, aln = load_alnfile(alnfile, reffile)
    alldepth = defaultdict(int)
    totallen = 0
    with open(f'{outprefix}.stats', 'w') as f:
        f.write('contig\tstart\tend\t%s\n' % ('\t'.join([str(x)+'X' for x in cutoffs])))
        for contig, start, stop in regions:
            print(contig, start, stop)
            contiglen = stop - start
            totallen += contiglen
            depths = count_cov(aln, contig, start, stop)
            nsites = []
            for cutoff in cutoffs:
                nsite = np.sum(depths >= cutoff)
                alldepth[cutoff] += nsite
                nsites.append(nsite)
            f.write('%s\t%s\t%s\t%s\n' % (contig, start, stop, '\t'.join([str(x) for x in nsites])))
    with open(f'{outprefix}.total', 'w') as f:
        f.write('totallen\t%s\n' % ('\t'.join([str(x)+'X' for x in cutoffs])))
        f.write('%s\t%s\n' % (totallen, '\t'.join([str(alldepth[cutoff]) for cutoff in cutoffs])))


if __name__ == '__main__':
    main()