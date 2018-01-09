# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 22:49:14 2018

@author: Caiyd
"""

import click
import numpy as np
import pandas as pd
from pysam import VariantFile



def ann2df(varin, sites: list):
    # 各个subfileds的解释http://snpeff.sourceforge.net/SnpEff_manual.html#input
    ann_subfields = ['Allele', 'Annotation', 'Annotation_Impact',
                      'Gene_Name', 'Gene_ID', 'Feature_Type', 'Feature_ID',
                      'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p', 'cDNA.pos/cDNA.length',
                      'CDS.pos/CDS.length', 'AA.pos/AA.length', 'Distance', 'ERRORS/WARNINGS/INFO']
    ann_records = []
    for contig, loc in sites:
        ann_record = []
        try:
            rec = next(varin.fetch(contig=contig, start=loc-1, stop=loc))
            for subfiled in zip(*[x.split('|') for x in rec.info['ANN']]):
                ann_record.append('|'.join(subfiled))
        except StopIteration: # 如果query位点不包含于注释信息中，调用next会引发StopIteration异常
            ann_record = [''] * 16
        ann_records.append(ann_record)
    return pd.DataFrame(ann_records, columns=ann_subfields)

@click.command()
@click.option('--regionfile', help='待注释的包含snp位置信息的文件，接受压缩格式, 第一行为header')
@click.option('--chr-col', help='序列名的列名')
@click.option('--chrom', help='使用固定的染色体号(regionfile只包含一条染色体且不含序列名的列)', type=str, default=None)
@click.option('--loc-col', help='snp坐标位置的列名, 单点的')
@click.option('--start-col', help='区间的起始位置, 与--loc-col冲突')
@click.option('--end-col', help='区间的终止位置, 与--loc-col冲突')
@click.option('--efffile', help='snpeff的注释结果,vcf(.gz)或bcf')
@click.option('--outprefix', help='输出后的注释后文件名前缀')
@click.option('--chunksize', help='每次处理nchunk行数据,默认10000行', default=10000, type=int)
def main(regionfile, chr_col, chrom, loc_col, start_col, end_col, efffile, outprefix, chunksize):
    varin = VariantFile(efffile)
    reader = pd.read_csv(regionfile, sep='\s+', iterator=True)
    loop = True
    nchunk = 1
    while loop:
        try:
            chunk = reader.get_chunk(chunksize).reset_index() # reset index是为了方便和后面注释的dfconcat
            print(f'{chunksize * nchunk} lines')
            if chrom:
                chunk['chrom'] = chrom
                sites = chunk[['chrom', loc_col]].values
            else:
                chunk[chr_col] = chunk[chr_col].astype(str)
                sites = chunk[[chr_col, loc_col]].values
            anndf = ann2df(varin, sites)
            outdf = pd.concat([chunk, anndf], axis=1).drop(columns='index')
            if nchunk != 1:
                outdf.to_csv(f'{outprefix}.ann.gz', sep='\t', index=False, header=False, compression='gzip', mode='ab')
            else:
                outdf.to_csv(f'{outprefix}.ann.gz', sep='\t', index=False, compression='gzip')
            nchunk += 1
        except StopIteration:
            loop = False


if __name__ == '__main__':
    main()

