# -*- coding: utf-8 -*-
"""
Created on Fri May 25 10:49:32 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import gzip
import click


def loaddepth(statfile, cutoff):
    sites = []
    with gzip.open(statfile, 'rb') as f:
        for line in f:
            tline = line.decode().strip().split()
            if float(tline[3]) >= cutoff:
                sites.append(f'{tline[0]}:{tline[2]}')
    return set(sites)



@click.command()
@click.option('--statfile', help='mosdepth的输出结果(regions.bed.gz)')
@click.option('--vcffile', help='vcf或bcf文件')
@click.option('--indiv', help='目标个体')
@click.option('--cutoff', help='大于等于cutoff的位点被保留', type=int)
@click.option('--outprefix', help='输出前缀')
def main(statfile, vcffile, indiv, cutoff, outprefix):
    """
    提取在mosdepth统计结果<statfile>中深度大于<cutoff>且在vcf文件<vcffile>中改个体<indiv>被成功分型(不为./.)的位点
    """
    sites = loaddepth(statfile, cutoff)
    with gzip.open(f'{outprefix}.{cutoff}x.sites.gz', 'wb') as f1:
        with gzip.open(vcffile, 'rb') as f2:
            for line in f2:
                tline = line.decode().strip().split('\t')
                if tline[0][0] != '#':
                    if f'{tline[0]}:{tline[1]}' in sites:
                        if tline[indiv_index].split(':')[0] != './.':
                            outline = f'{tline[0]}\t{tline[1]}\n'
                            f1.write(outline.encode())
                elif tline[0] == '#CHROM':
                    indiv_index = tline.index(indiv)

if __name__ == '__main__':
    main()
