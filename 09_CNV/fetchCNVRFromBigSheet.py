# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 19:26:46 2017

@author: Caiyd
"""

import sys
import pysam
import click
import numpy as np


def fetchcnvr(tbx, region):
    cnvay = []
    try:
        for row in tbx.fetch(region):
            cnvay.append(row.strip().split()[3:])
    except ValueError:
        pass
    cnvay = np.array(cnvay, dtype=float)
    if len(cnvay)>0:
        return np.median(cnvay, axis=0).tolist() # 取中位数
    else:
        return None


def loadregion(regionfile):
    regionlist = [x.strip() for x in open(regionfile).readlines()]
    return regionlist



@click.command()
@click.argument('bigsheet')
@click.argument('regionfile')
@click.argument('outfile')
def main(bigsheet, regionfile, outfile):
    if float(pysam.__version__) < 0.13:
        sys.exit('pysam version need >= 0.13')
    tbx = pysam.TabixFile(bigsheet)
    smlist = tbx.header[0].split('\t')[3:]
    regionlist = loadregion(regionfile)
    with open(outfile, 'w') as f:
        f.write('region\t' + '\t'.join(smlist) + '\n')
        with click.progressbar(regionlist) as bar:
            for region in bar:
                cnvlist = fetchcnvr(tbx, region)
                if cnvlist:
                    f.write(region + '\t' + '\t'.join([str(x) for x in cnvlist]) + '\n')


if __name__ == '__main__':
    main()

