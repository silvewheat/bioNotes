# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 21:09:53 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import sys
import click
import pysam
import pandas as pd
from collections import Counter


def load_alnfile(alnfile, reffile):
    if alnfile.split('.')[-1] == 'cram':
        if not reffile:
            sys.exit("cram file need refrence fasta in --reffile")
        aln = pysam.AlignmentFile(alnfile, reference_filename=reffile)
    else:
        aln = pysam.AlignmentFile(alnfile)
    smid = aln.header['RG'][0]['SM']
    return smid, aln


def get_clip_site(record, cigarstring):
    """
    如果存在多个clip位点，只报告处于最右边或者最左边的clip位点（右>左）
    """
    cigarstring = cigarstring.strip('0123456789')
    blocks = record.get_blocks()
    if (cigarstring[-1] == 'S') or (cigarstring[-1] == 'H'):
        edge_left = blocks[-1][-1] # covered
        edge_right = edge_left + 1 # uncovered
        return 'rightclip', edge_left, edge_right
    elif (cigarstring[0] == 'S') or (cigarstring[0] == 'H'):
        edge_right = blocks[0][0] # covered
        edge_left = edge_right - 1 # uncovered
        return 'leftclip', edge_left, edge_right
    else:
        return 'MultiClip', None, None





@click.command()
@click.option('--alnfile', help='cram/bam/sam的文件')
@click.option('--reffile', help='输入的参考基因文件,cram文件必须要,bam不用', default=None)
@click.option('--region', help='regionstring, like 1:100-200')
@click.option('--outfile', help='output file')
def main(alnfile, reffile, region, outfile):
    """
    \b
    统计cover到--region的reads和clip位点
    leftclip: reads covered edge_right, uncovered edge_left
    righrclip: reads covered edge_left, uncovered edge_right
    """
    result = []
    smid, aln = load_alnfile(alnfile, reffile)
    for record in aln.fetch(region=region):
        cigarstring = record.cigarstring
        count = Counter(cigarstring)
        n_soft = count.get('S', 0)
        n_hard = count.get('H', 0)
        if n_soft + n_hard == 0:
            continue
        else:
            clip_type, edge_left, edge_right = get_clip_site(record, cigarstring)
            name = record.query_name
            result.append([smid, name, clip_type, edge_left, edge_right])
    df = pd.DataFrame(result, columns=['sample_id', 'read_name', 'clip_type', 'edge_left', 'edge_right'])
    df.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()




