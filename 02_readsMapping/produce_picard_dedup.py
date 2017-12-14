# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 10:31:20 2017

@author: Caiyd
"""


import os
import click
import pysam
from collections import defaultdict


def produce_cmd(ref, groupbams, outmetrics, outcram, picard, tmpdir, mem, samtools):
    cmd = f"""java -Djava.io.tmpdir={tmpdir} \\
    -Xmx{mem}g \\
    -jar {picard} \\
    MarkDuplicates \\
"""
    nbams = len(groupbams)
    for i, inbam in enumerate(groupbams, 1):
        if i != nbams:
            cmd += f"""    I={inbam} \\\n"""
        else:
            cmd += f"""    I={inbam} \\"""
    cmd += f"""
    O=/dev/stdout \\
    M={outmetrics} \\
    REMOVE_DUPLICATES=true \\
    VALIDATION_STRINGENCY=LENIENT \\
    COMPRESSION_LEVEL=0 | \\
{samtools} \\
    view \\
    -T {ref} \\
    -C \\
    -o {outcram} \\
    -
    """
    return cmd


@click.command()
@click.option('--ref', help='参考基因组的路径')
@click.option('--outdir', help='输出去冗余结果的路径')
@click.option('--picard', help='picard的路径')
@click.option('--tmpdir', help='输出临时文件的目录')
@click.option('--mem', help='运行picard的内存，单位是GB')
@click.option('--samtools', help='samtools的路径', default='samtools')
@click.argument('inbams', nargs=-1)
def main(ref, outdir, picard, tmpdir, mem, samtools, inbams):
    groups = defaultdict(list)
    for inbam in inbams:
        f = pysam.AlignmentFile(inbam)
        smid = f.header['RG'][0]['SM']
        groups[smid].append(inbam)
    for ngroup, (smid, groupbams) in enumerate(groups.items(), 1):
        print(ngroup, smid, groupbams)
        outcram = os.path.join(outdir, f'{smid}.cram')
        outmetrics = os.path.join(outdir, f'{smid}.marked_dup_metrics.txt')
        cmd = produce_cmd(ref, groupbams, outmetrics, outcram, picard, tmpdir, mem, samtools)
        with open(f'dedup_{ngroup}.sh', 'w') as f:
            f.write(cmd)


if __name__ == '__main__':
    main()



