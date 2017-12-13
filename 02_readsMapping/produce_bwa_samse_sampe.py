# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 21:28:46 2017

@author: Caiyd
"""

import os
import click


def produce_samse(ref, infastq, insai, bwa, samtools, outfile, rgid, lb, pl, smid, outfmt):
    return f"""{bwa} \\
    samse \\
    -r '@RG\\tID:{rgid}\\tLB:{lb}\\tPL:{pl}\\tSM:{smid}' \\
    {ref} \\
    {insai} \\
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
echo "samse is done"
    """

def produce_sampe(ref, in1fastq, in2fastq, in1sai, in2sai, bwa, samtools, outfile, rgid, lb, pl, smid, outfmt):
    return f"""{bwa} \\
    sampe \\
    -r '@RG\\tID:{rgid}\\tLB:{lb}\\tPL:{pl}\\tSM:{smid}' \\
    {ref} \\
    {in1sai} \\
    {in2sai} \\
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



def split_sepe(fastqs):
    fastqs = list(fastqs)
    fastqs.sort()
    sefiles = []
    pefiles = []
    for fastq in fastqs:
        basename, suffix = fastq.split('.', 1)
        if basename[-2:] == '_1':
            if (basename[:-2] + '_2.' + suffix) in fastqs: # 如果有对应_1的_2
                pefiles.append(fastq)
                pefiles.append(basename[:-2] + '_2.' + suffix)
        elif basename[-2:] == '_2':
            if fastq not in pefiles: # fastqs拍过序，如果是成对的，_1一定在_2之前， 这会儿还没存就是单端了
                sefiles.append(fastq)
        else:
            sefiles.append(fastq)
    return sefiles, pefiles


@click.command()
@click.option('--ref', help='参考基因组的路径')
@click.option('--outdir', help='输出比对结果的文件')
@click.option('--insaidir', help='sai文件所在的文件夹')
@click.option('--bwa', help='BWA的路径', default='bwa')
@click.option('--samtools', help='samtools的路径', default='samtools')
@click.option('--outfmt', help='输出文件的格式,BAM或CRAM,默认BAM', default='BAM')
@click.option('--pl', help='测序平台, default is ILLUMINA', default='ILLUMINA')
@click.argument('infastqs', nargs=-1)
def main(ref, outdir, insaidir, bwa, samtools, outfmt, pl, infastqs):
    sefiles, pefiles = split_sepe(infastqs)
    for nfile, infastq in enumerate(sefiles, 1):
        basename = os.path.basename(infastq).split('.')[0]
        insai = os.path.join(insaidir, basename + '.sai')
        outfile = os.path.join(outdir, f'{basename}.{outfmt.lower()}')
        rgid = basename
        lb = basename
        smid = basename.split('_')[0]
        cmd = produce_samse(ref, infastq, insai, bwa, samtools, outfile, rgid, lb, pl, smid, outfmt)
        with open(f'bwa_samse_{nfile}.sh', 'w') as f:
            f.write(cmd)
    for npair, (in1fastq, in2fastq) in enumerate(zip(pefiles[::2], pefiles[1::2]), 1):
        basename1 = os.path.basename(in1fastq).split('.')[0]
        basename2 = os.path.basename(in2fastq).split('.')[0]
        in1sai = os.path.join(insaidir, basename1 + '.sai')
        in2sai = os.path.join(insaidir, basename2 + '.sai')
        outfile = os.path.join(outdir, f'{basename1[:-2]}.{outfmt.lower()}')
        rgid = basename1[:-2]
        lb = basename1[:-2]
        smid = basename1.split('_')[0]
        cmd = produce_sampe(ref, in1fastq, in2fastq, in1sai, in2sai, bwa, samtools, outfile, rgid, lb, pl, smid, outfmt)
        with open(f'bwa_sampe_{npair}.sh', 'w') as f:
            f.write(cmd)



if __name__ == '__main__':
    main()
