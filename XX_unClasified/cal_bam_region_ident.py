# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 15:51:57 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import sys
import click
import pysam
import numpy as np
import pandas as pd


def load_region(regionfile):
    regions = []
    with open(regionfile) as f:
        for line in f:
            if line[0] != '#':
                sline = line.strip().split()
                regions.append([sline[0], int(sline[1]), int(sline[2])])
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


def cal_record_ident(region, aln):
    idents = []
    contig, start, stop = region
    for record in aln.fetch(contig=contig, start=start, stop=stop):
        ident = 1 - (record.get_tag('NM') / record.query_alignment_length)
        idents.append(ident)
    if idents:
        return np.mean(idents), np.median(idents), np.std(idents), np.min(idents), np.max(idents), len(idents)
    else:
        return None, None, None, None, None, None




@click.command()
@click.option('--alnfile', help='cram/bam/sam的文件')
@click.option('--reffile', help='输入的参考基因文件,cram文件必须要,bam不用', default=None)
@click.option('--regionfile', help='需要计算的区域，bed格式, 如有header,开头为#')
@click.option('--outfile', help='输出结果文件的前缀')
def main(alnfile, reffile, regionfile, outfile):
    """
    计算比对到--regionfile中
    """
    result = []
    smid, aln = load_alnfile(alnfile, reffile)
    regions = load_region(regionfile)
    for region in regions:
        mean, median, std, imin, imax, nreads = cal_record_ident(region, aln)
        result.append(region + [smid, mean, median, std, imin, imax, nreads])
    df = pd.DataFrame(result,
                      columns=['chrom', 'start', 'end', 'sample_id', 'ident_mean', 'ident_median', 'ident_std',
                               'ident_min', 'ident_max', 'num_reads'])
    df.to_csv(outfile, sep='\t', index=False, na_rep='NaN')


if __name__ == '__main__':
    main()


