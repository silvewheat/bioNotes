# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 15:53:28 2017

@author: Caiyd
"""

import os
import click



def produce_RTC(ref, gatk, mem, inbam, outinterval, tmpdir, outbam, nt):
    return f"""java -Djava.io.tmpdir={tmpdir} \\
    -Xmx{mem * nt}g \\
    -jar {gatk} \\
    -T RealignerTargetCreator \\
    -R {ref} \\
    -I {inbam} \\
    -nt {nt} \\
    -o {outinterval} \\
    -allowPotentiallyMisencodedQuals \\
    && \\
    echo 'RealignerTargetCreator Done'
"""

def produce_IR(ref, gatk, mem, inbam, outinterval, tmpdir, outbam, samtools):
    return f"""java -Djava.io.tmpdir={tmpdir} \\
    -Xmx{mem}g \\
    -jar {gatk} \\
    -T IndelRealigner \\
    -R {ref} \\
    -I {inbam} \\
    -targetIntervals {outinterval} \\
    -allowPotentiallyMisencodedQuals \\
    -o {outbam} \\
    && \\
    echo 'IndelRealigner Done' \\
    && \\
{samtools} \\
    index \\
    {outbam}
"""

@click.command()
@click.option('--ref', help='参考基因组的路径')
@click.option('--outdir', help='输出局部重比后的cram文件路径')
@click.option('--gatk', help='gatk的路径')
@click.option('--samtools', help='samtools的路径')
@click.option('--tmpdir', help='输出临时文件的目录')
@click.option('--mem', type=int, help='单线程的java启动内存，单位GB，默认4', default=4)
@click.option('--nt', type=int, help='运行RTC时的线程数,默认4', default=4)
@click.argument('inbams', nargs=-1)
def main(ref, outdir, gatk, samtools, tmpdir, mem, nt, inbams):
    for n, inbam in enumerate(inbams, 1):
        smid = os.path.basename(inbam).split('.')[0]
        outbam = os.path.join(outdir, f'{smid}.realn.bam')
        outinterval = os.path.join(outdir, f'{smid}.RTC.intervals')
        RTC_cmd = produce_RTC(ref, gatk, mem, inbam, outinterval, tmpdir, outbam, nt)
        IR_cmd = produce_IR(ref, gatk, mem, inbam, outinterval, tmpdir, outbam, samtools)
        with open(f'RTC_{n}.sh', 'w') as f:
            f.write(RTC_cmd)
        with open(f'IR_{n}.sh', 'w') as f:
            f.write(IR_cmd)


if __name__ == '__main__':
    main()
