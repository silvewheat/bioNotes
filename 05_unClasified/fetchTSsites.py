# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 11:18:22 2018

@author: Caiyd
"""

import click
import gzip
from pysam import VariantFile


@click.command()
@click.option('--varfile', help='vcf(.gz), bcf文件')
@click.option('--outprefix', help='{outprefix}.TS.sites.gz')
def main(varfile, outprefix):
    """
    从vcf(.gz), bcf文件中提取转换（transition）的位点
    保存在{outprefix}.TS.sites.gz文件中
    """
    transitions = {'AT', 'TA', 'CG', 'GC'}
    varin = VariantFile(varfile)
    with gzip.open(f'{outprefix}.TS.sites.gz', 'wb') as f:
        for rec in varin.fetch():
            ref = rec.ref
            alts = rec.alts
            if len(alts) == 1:
                alt = alts[0]
            else:
                continue
            if f'{ref}{alt}' in transitions:
                outstr = f'{rec.chrom}\t{rec.pos}\n'.encode()
                f.write(outstr)


if __name__ == '__main__':
    main()
