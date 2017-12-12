# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 11:12:38 2017

@author: Caiyd
"""

import os
import click


def produce_cmd(ref, infastq, outsai, bwa, nt):
    return f"""{bwa} \\
    aln \\
    -t {nt} \\
    -l 1024 \\
    {ref} \\
    {infastq} \\
    > {outsai} \\
    && \\
    echo "bwa aln is done!"
    """


@click.command()
@click.option('--ref', help='参考基因组的路径')
@click.option('--outsaidir', help='输出sai结果的文件夹')
@click.option('--nt', help='线程数', default=1)
@click.option('--bwa', help='BWA的路径', default='bwa')
@click.argument('infastqs', nargs=-1)
def main(ref, outsaidir, nt, bwa, infastqs):
    for nfile, infastq in enumerate(infastqs, 1):
        basename = os.path.basename(infastq).split('.')[0]
        saibasename = basename + '.sai'
        outsai = os.path.join(outsaidir, saibasename)
        cmd = produce_cmd(ref, infastq, outsai, bwa, nt)
        with open(f'bwa_aln_{nfile}.sh', 'w') as f:
            f.write(cmd)


if __name__ == '__main__':
    main()

