# -*- coding: utf-8 -*-
"""
Created on Sat Apr 14 20:47:39 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np
import pandas as pd
import h5py
import allel
from itertools import combinations


def getAC(genotypes, variant_selection, sample_selection):
    return genotypes.subset(variant_selection, sample_selection).count_alleles()


def cal_all_ac(genotypes, samples, allsamples, variant_selection):
    """
    cal all single sample allel count
    """
    print('cal single sample allele count')
    nsample = len(samples)
    ac_dict = {}
    n = 0
    for sample in samples:
        n += 1
        print(f'{n}/{nsample}\t{sample}')
        sample_selection =  [True if x == sample else False for x in allsamples]
        ac_dict[sample] = getAC(genotypes, variant_selection, sample_selection)
    print('cal single sample allele count finish')
    return ac_dict


def cal_outgroup_ac(genotypes, outgroupfile, allsamples, variant_selection):
    """
    cal outgroup ac
    """
    print('cal outgroup allel count')
    outgroup = [x.strip() for x in open(outgroupfile).readlines()]
    outgroup_selection = [True if x in outgroup else False for x in allsamples]
    ac_outgroup = getAC(genotypes, variant_selection, outgroup_selection)
    print('cal outgroup allel count done')
    return ac_outgroup


@click.command()
@click.option('--h5file', help='h5file, produced from vcf')
@click.option('--samplesfile', help='samples id used to cal pairwise outgroup-f3')
@click.option('--outgroupfile', help='outgroup sample id file')
@click.option('--outprefix', help='outfile prefix')
@click.option('--blen', help='Block size (number of variants use in block-jackknife), default is 10_000', default=10000)
def main(h5file, samplesfile, outgroupfile, outprefix, blen):
    """
    use samples in samplesfile to calculate pair-wise outgroup-f3
    use samples in outgroupfile as outgroup
    h5file generate from vcffile from scikit-allele(1.1.10)
    import allel; allel.vcf_to_hdf5('in.vcf.gz', 'out.h5')
    """
    print(__doc__)
    print('scikit-allel', allel.__version__)
    samples = [x.strip() for x in open(samplesfile).readlines()] # 待计算个体
    callset = h5py.File(h5file, mode='r')
    allsamples = list(callset['samples']) # vcf包含的全部个体
    calldata = callset['calldata']
    genotypes = allel.GenotypeChunkedArray(calldata['GT'])
    variant_selection = np.full((genotypes.shape[0]+1), True) # 选择vcf中的全部位点
    ac_outgroup = cal_outgroup_ac(genotypes, outgroupfile, allsamples, variant_selection)
    ac_dict = cal_all_ac(genotypes, samples, allsamples, variant_selection)
    print('begin to cal outgroup f3')
    n_comb = len(list(combinations(samples, 2)))
    print(f'total combinations is {n_comb}')
    n_iter = 0
    n_samples = len(samples)
    f3_ay = np.full((n_samples, n_samples), None)
    z_ay = np.full((n_samples, n_samples), None)
    for sample1, sample2 in combinations(samples, 2):
        x = samples.index(sample1)
        y = samples.index(sample2)
        n_iter += 1
        print(f'{n_iter}/{n_comb}')
        f3, se, z, vb, vj = allel.average_patterson_f3(ac_dict[sample1], ac_dict[sample2], ac_outgroup, blen)
        f3_ay[x, y] = f3
        f3_ay[y, x] = f3
        z_ay[x, y] = z
        z_ay[y, x] = z
    pd.DataFrame(f3_ay, columns=samples, index=samples).to_csv(f'{outprefix}.f3.tsv', sep='\t')
    pd.DataFrame(z_ay, columns=samples, index=samples).to_csv(f'{outprefix}.z.tsv', sep='\t')



if __name__ == '__main__':
    main()







