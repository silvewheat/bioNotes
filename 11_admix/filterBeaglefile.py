# -*- coding: utf-8 -*-
"""
Created on Wed May 16 10:53:26 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import gzip
import click


@click.command()
@click.option('--sites', help='一列, gz压缩, chrom:pos')
@click.option('--gzbeaglefile', help='gz压缩beagle文件')
@click.option('--outfile', help='输出文件，后缀是.gz')
def main(sites, gzbeaglefile, outfile):
    sites = set([x.decode().strip() for x in gzip.open(sites, 'rb').readlines()])
    print(f'load {len(sites)} sites')
    with gzip.open(outfile, 'wb') as f1:
        with gzip.open(gzbeaglefile, 'rb') as f2:
            header = f2.readline()
            f1.write(header)
            for line in f2:
                site = line.decode().split()[0]
                if site in sites:
                    f1.write(line)

if __name__ == '__main__':
    main()

