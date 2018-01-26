# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 21:02:34 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import sys
import click
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pysam



def get_ref_len(alnfile, contig):
    for n, seqinfo in enumerate(alnfile.header['SQ']):
        if seqinfo['SN'] == contig:
            return seqinfo['LN']


def get_contig_mapped(alnfile, contig):
    """
    返回一条contig上mapping上的reads综述
    """
    for record in alnfile.get_index_statistics():
        if record.contig == contig:
            count = record.mapped
            if count > 0: # cram 文件的 index 中似乎不包含比对上的 reads 数
                return count
            else:
                count = alnfile.count(contig)
                return count
    else:
        sys.exit(f"no contig named {contig}")


def get_median(file, contig, reffile, winsize=1000000, stepsize=1000000):
    """
    返回个窗口内的reads count中位数
    """
    if reffile:
        alnfile = pysam.AlignmentFile(file, reference_filename=reffile)
    else:
        alnfile = pysam.AlignmentFile(file)
    smid = alnfile.header['RG'][0]['SM']
    contiglen = get_ref_len(alnfile, contig)
    counts = []
    for start in range(1, contiglen, stepsize):
        counts.append(alnfile.count(contig, start, start+winsize-1))
    return smid, np.median(counts)


def draw(ratios, outprefix):
    ratios = np.array(ratios)
    ratios[ratios > 3] = 3
    fig, ax = plt.subplots(1, 1, figsize=(8,4))
    bins = ax.hist(ratios[ratios<=3], bins = 50)
    ax.set_xlim(0.5, 2.5)
    ax.vlines(1.5, 0, max([y for x in bins[:-1] for y in x]), linestyles='dashed')
    ax.set_xticks(np.arange(0.5, 2.6, 0.5))
    ax.set_xticklabels(['≤0.5', '1.0', '1.5', '2.0', '≥2.5'])
    plt.savefig(f'{outprefix}.pdf')
    plt.close()


@click.command()
@click.option('--filelist', help='cram/bam/sam的文件列表,一行一个文件,带路径')
@click.option('--reffile', help='输入的参考基因文件,cram文件必须要,bam不用', default=None)
@click.option('--chr1-id', help='1号染色体在参考基因中的名字')
@click.option('--chrx-id', help='x染色体在参考基因组中的名字')
@click.option('--outprefix', help='输出结果文件的前缀')
def main(filelist, reffile, chr1_id, chrx_id, outprefix):
    """
    \b
    根据CRAM/BAM文件中常染色体和X染色体的reads count比例来判断个体性别
    输出的chr1_count和chrX_count分别是chr1和chrX上的总reads数
    chr1/chrX是前面的reads数用各自染色体长度标准化后相除的比值
    1.5作为雌雄的分界线
    """
    files = [x.strip() for x in open(filelist).readlines()]
    ratios = []
    with open(f'{outprefix}.stats', 'w') as f:
        f.write('smid\tchr1_count\tchrX_count\tchr1_scaled\tchrX_scaled\tchr1/chrX\tsex\n')
        for file in files:
            if file.split('.')[-1] == 'cram':
                if not reffile:
                    sys.exit("cram file need refrence fasta in --reffile")
                alnfile = pysam.AlignmentFile(file, reference_filename=reffile)
            else:
                alnfile = pysam.AlignmentFile(file)
            smid = alnfile.header['RG'][0]['SM']
            chr1_count = get_contig_mapped(alnfile, chr1_id)
            chr1_len = get_ref_len(alnfile, chr1_id)
            chrx_count = get_contig_mapped(alnfile, chrx_id)
            chrx_len = get_ref_len(alnfile, chrx_id)
            if chrx_count > 0:
                ratio = (chr1_count/chr1_len) / (chrx_count/chrx_len)
                ratios.append(ratio)
                sex = 'male' if ratio >= 1.5 else 'female'
            else:
                ratio = 'nan'
                sex = 'unknow'
            chr1_scaled = int((chr1_count/chr1_len)*((chr1_len+chrx_len)/2))
            chrx_scaled = int((chrx_count/chrx_len)*((chr1_len+chrx_len)/2))
            f.write(f'{smid}\t{chr1_count}\t{chrx_count}\t{chr1_scaled}\t{chrx_scaled}\t{sex}\n')
    draw(ratios, outprefix)


if __name__ == '__main__':
    main()
