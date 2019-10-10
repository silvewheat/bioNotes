# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 11:12:38 2017
@author: Caiyd
"""

import os
import click


def produce_mempe(ref, in1fastq, in2fastq, bwa, samtools, outfile, rgid, lb, pl, smid, outfmt, nt):
    return f"""{bwa} \\
    mem \\
    -t {nt} \\
    -R '@RG\\tID:{rgid}\\tLB:{lb}\\tPL:{pl}\\tSM:{smid}' \\
    {ref} \\
    {in1fastq} \\
    {in2fastq} | \\
{samtools} \\
    fixmate \\
    -r \\
    - - | \\
{samtools} \\
    sort \\
    --reference {ref} \\
    --output-fmt {outfmt} \\
    -o {outfile} \\
    - \\
    && \\
echo "sampe is done"
    """

def produce_memse(ref, infastq, bwa, samtools, outfile, rgid, lb, pl, smid, outfmt, nt):
    return f"""{bwa} \\
    mem \\
    -t {nt} \\
    -R '@RG\\tID:{rgid}\\tLB:{lb}\\tPL:{pl}\\tSM:{smid}' \\
    {ref} \\
    {infastq} | \\
{samtools} \\
    fixmate \\
    -r \\
    - - | \\
{samtools} \\
    sort \\
    --reference {ref} \\
    --output-fmt {outfmt} \\
    -o {outfile} \\
    - \\
    && \\
echo "sampe is done"
    """


def split_sepe(fastqs):
    fastqs = list(fastqs)
    fastqs.sort()
    print(fastqs)
    sefiles = []
    pefiles = []
    for fastq in fastqs:
        path, file = os.path.split(fastq)
        basename, suffix = file.split('.', 1)
        print(basename, suffix)
        print(basename[:-2])
        if basename[-2:] == '_1':
            if (os.path.join(path, basename[:-2] + '_2.' + suffix)) in fastqs: # 如果有对应_1的_2
                pefiles.append(fastq)
                pefiles.append(os.path.join(path, basename[:-2] + '_2.' + suffix))
        elif basename[-2:] == '_2':
            if fastq not in pefiles: # fastqs拍过序，如果是成对的，_1一定在_2之前， 这会儿还没存就是单端了
                sefiles.append(fastq)
        else:
            sefiles.append(fastq)
    return sefiles, pefiles


@click.command()
@click.option('--ref', help='参考基因组的路径')
@click.option('--outdir', help='输出比对结果的文件')
@click.option('--bwa', help='BWA的路径', default='bwa')
@click.option('--samtools', help='samtools的路径', default='samtools')
@click.option('--outfmt', help='输出文件的格式,BAM或CRAM,默认BAM', default='BAM')
@click.option('--pl', help='测序平台, default is ILLUMINA', default='ILLUMINA')
@click.option('--nt', help='线程数', default=4)
@click.argument('infastqs', nargs=-1)
def main(ref, outdir, bwa, samtools, outfmt, pl, nt, infastqs):
    "python produce_bwa_mem.py [options] /home/data/*.fq.gz"
    sefiles, pefiles = split_sepe(infastqs)
    for nfile, infastq in enumerate(sefiles, 1):
        basename = os.path.basename(infastq).split('.')[0]
        outfile = os.path.join(outdir, f'{basename}.{outfmt.lower()}')
        rgid = basename
        lb = basename
        smid = basename.split('_')[0]
        cmd = produce_memse(ref, infastq, bwa, samtools, outfile, rgid, lb, pl, smid, outfmt, nt)
        with open(f'bwa_memse_{nfile}.sh', 'w') as f:
            f.write(cmd)
    for npair, (in1fastq, in2fastq) in enumerate(zip(pefiles[::2], pefiles[1::2]), 1):
        basename1 = os.path.basename(in1fastq).split('.')[0]
        basename2 = os.path.basename(in2fastq).split('.')[0]
        assert basename1[:-2] == basename2[:-2]
        outfile = os.path.join(outdir, f'{basename1[:-2]}.{outfmt.lower()}')
        rgid = basename1[:-2]
        lb = basename1[:-2]
        smid = basename1.split('_')[0]
        cmd = produce_mempe(ref, in1fastq, in2fastq, bwa, samtools, outfile, rgid, lb, pl, smid, outfmt, nt)
        with open(f'bwa_mempe_{npair}.sh', 'w') as f:
            f.write(cmd)



if __name__ == '__main__':
    main()
