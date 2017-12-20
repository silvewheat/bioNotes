# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 15:53:28 2017

@author: Caiyd
"""

import os
import click



def produce_RTC(ref, gatk, mem, inbam, outinterval, tmpdir, outcram, nt):
    return f"""java -Djava.io.tmpdir={tmpdir} \\
    -Xmx{mem}g \\
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

def produce_IR(ref, gatk, mem, inbam, outinterval, tmpdir, outcram):
    return f"""java -Djava.io.tmpdir={tmpdir} \\
    -Xmx{mem}g \\
    -jar {gatk} \\
    -T IndelRealigner \\
    -R {ref} \\
    -I {inbam} \\
    -targetIntervals {outinterval} \\
    -o {outcram} \\
    && \\
    echo 'IndelRealigner Done'
"""

@click.command()
@click.option('--ref', help='参考基因组的路径')
@click.option('--outdir', help='输出局部重比后的cram文件路径')
@click.option('--gatk', help='gatk的路径')
@click.option('--tmpdir', help='输出临时文件的目录')
@click.option('--mem', help='运行gatk的内存，单位是GB')
@click.option('--nt', help='运行RTC时的线程数')
@click.argument('inbams', nargs=-1)
def main(ref, outdir, gatk, tmpdir, mem, nt, inbams):
    for n, inbam in enumerate(inbams, 1):
        smid = os.path.basename(inbam).split('.')[0]
        outcram = os.path.join(outdir, f'{smid}.cram')
        outinterval = os.path.join(outdir, f'{smid}.RTC.intervals')
        RTC_cmd = produce_RTC(ref, gatk, mem, inbam, outinterval, tmpdir, outcram, nt)
        IR_cmd = produce_IR(ref, gatk, mem, inbam, outinterval, tmpdir, outcram)
        with open(f'RTC_{n}.sh', 'w') as f:
            f.write(RTC_cmd)
        with open(f'IR_{n}.sh', 'w') as f:
            f.write(IR_cmd)


if __name__ == '__main__':
    main()
