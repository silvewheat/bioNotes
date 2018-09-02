# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 16:09:30 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import sys
import click
import numpy as np
import pysam




def print_flush(text):
    sys.stdout.write(' ' * 30 + '\r')
    sys.stdout.flush()
    sys.stdout.write(f' > {text}' + '\r')
    sys.stdout.flush()


def produce_record(varfile, samples, nsites):
    var = pysam.VariantFile(varfile)
    var.subset_samples(samples)
    for nrec, record in enumerate(var.fetch(), 1):
        print_flush(f'{nrec} / {nsites}')
        gts = np.array([x['GT'] for x in record.samples.values()])
        n_ref = np.sum(gts == 0)
        n_alt = np.sum(gts == 1)
        n_hap = n_ref + n_alt
        yield f'{nrec}\t{n_hap}\t2\t{n_ref}\t{n_alt}\n'


@click.command()
@click.option('--varfile', help='vcf(.gz)/bcf file')
@click.option('--nsites', help='number of sites in varfile', type=int)
@click.option('--poplists', help='files contained sample IDs in each pop', multiple=True)
@click.option('--outfile', help='outfile name')
def main(varfile, nsites, poplists, outfile):
    """
    \b
    convert vcf(.gz)/bcf (biallelic only) to BayeScan (Codominant markers format)
    http://cmpg.unibe.ch/software/BayeScan/index.html
    e.g.
    python vcf2BayeScan.py --varfile in.vcf.gz --nsites 10000 --poplists pop1.id --poplists pop2.id --outfile out.txt
    poplists file contain only one column, one id per row.
    """
    npop = len(poplists)
    with open(outfile, 'w') as f:
        f.write(f'[loci]={nsites}\n\n[populations]={npop}\n\n')
        for npop, popfile in enumerate(poplists, 1):
            print(f'pop{npop}: {popfile}')
            f.write(f'[pop]={npop}\n')
            samples = [x.strip() for x in open(popfile)]
            for line in produce_record(varfile, samples, nsites):
                f.write(line)
            f.write('\n')


if __name__ == '__main__':
    main()





