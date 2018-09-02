# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 16:16:40 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import sys
import click
import numpy as np
import pysam


def load_alnfile(alnfile, reffile):
    if alnfile.split('.')[-1] == 'cram':
        if not reffile:
            sys.exit("cram file need refrence fasta in --reffile")
        aln = pysam.AlignmentFile(alnfile, reference_filename=reffile)
    else:
        aln = pysam.AlignmentFile(alnfile)
    smid = aln.header['RG'][0]['SM']
    return smid, aln


def count_cov(aln, region: str, quality: int):
    """
    region: region string
    返回一个列表，包含region中每个位点的深度
    """
    contig = ":".join(region.split(':')[:-1])
    pos = region.split(':')[-1]
    start, stop = [int(x) for x in pos.split('-')]
    return np.sum(aln.count_coverage(contig=contig, start=start, stop=stop, quality_threshold=quality), axis=0)


@click.command()
@click.option('--alnfile', help='BAM/CRAM/SAM')
@click.option('--reffile', help='refrence file (only for CRAM)', default=None)
@click.option('--queryregion', help='region need to count')
@click.option('--refregion', help='use depth in this region as refrence to normalize, 必需包含起始和终止')
@click.option('--quality', help='minimum quality score(in phred) a base has to reach to be counted, default is 0', default=0, type=int)
@click.option('--outfile', help='output file')
def main(alnfile, reffile, queryregion, refregion, quality, outfile):
    """
    \b
    统计alnfile中queryregion的reads深度，
    并使用refregion中的深度为参考进行标准化（refregion的深度作为1）
    """
    print(__doc__)
    smid, aln = load_alnfile(alnfile, reffile)
    refdps = count_cov(aln, refregion, quality)
    refdp = np.median(refdps)
    refmean = np.mean(refdps)
    print(f'refregion DP median: {refdp}, mean: {refmean}')
    querydp = count_cov(aln, queryregion, quality)
    querymean = np.mean(querydp)
    querymean_nm = querymean / refdp
    querymedian = np.median(querydp)
    querymedian_nm = querymedian / refdp
    with open(outfile, 'w')  as f:
        f.write('refDP\tquery_mean\tqueryDP_mean_NM\tquery_median\tquery_median_NM\n')
        f.write(f'{refdp}\t{querymean}\t{querymean_nm}\t{querymedian}\t{querymedian_nm}\n')

if __name__ == '__main__':
    main()

