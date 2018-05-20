# -*- coding: utf-8 -*-
"""
Created on Wed May 16 09:38:20 2018

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
@click.option('--cutoff', help='大于等于cutoff的位点被保留', type=int)
@click.option('--snpfile')
@click.option('--genofile')
@click.option('--outprefix')
def main(statfile, cutoff, snpfile, genofile, outprefix):
    sites = loaddepth(statfile, cutoff)
    print(f'load depth done, {len(sites)} sites remained')
    index = []
    with open(f'{outprefix}.snp', 'w') as f1:
        with open(snpfile) as f2:
            for nline, line in enumerate(f2):
                tline = line.strip().split()
                if tline[0] in sites:
                    f1.write(line)
                    index.append(nline)
    index = set(index)
    print(f'filter snpfile done, {len(index)} sites remained')
    with open(f'{outprefix}.geno', 'w') as f1:
        with open(genofile) as f2:
            for nline, line in enumerate(f2):
                if nline in index:
                    f1.write(line)

if __name__ == '__main__':
    main()

