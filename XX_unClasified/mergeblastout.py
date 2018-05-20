# -*- coding: utf-8 -*-
"""
Created on Thu May  3 09:33:00 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import os
import click



def loadblastout(blastout):
    with open(blastout) as f:
        records = []
        for line in f:
            if line[0] != '#':
                sline = line.strip().split()
                qid = sline[0]
                start = int(sline[6])
                end = int(sline[7])
                records.append(f'{qid}\t{min(start, end)}\t{max(start, end)}')
            else:
                if records:
                    tmpfile = f'tmp_{qid}.bed'
                    with open(tmpfile, 'w') as f_out:
                        f_out.write('\n'.join(records))
                        records = []
                    yield qid, calbedlen(tmpfile)
                    os.remove(tmpfile)

def calquerycov(blastout, querylen):
    queryaln = {}
    with open(blastout) as f:
        records = []
        for line in f:
            if line[0] != '#':
                sline = line.strip().split()
                qid = sline[0]
                start = int(sline[6])
                end = int(sline[7])
                records.append(f'{qid}\t{min(start, end)}\t{max(start, end)}')
            else:
                if records:
                    tmpfile = f'tmp_{qid}.bed'
                    with open(tmpfile, 'w') as f_out:
                        f_out.write('\n'.join(records))
                        records = []
                    queryaln[qid] = calbedlen(tmpfile) / querylen[qid]
                    os.remove(tmpfile)
    return queryaln


def calbedlen(bedfile):
    cmd = "sort -k2,2n %s | bedtools merge -d 1 | awk '{a+=($3-$2+1)};END{print a}'" % bedfile
    try:
        return int(os.popen(cmd).readline().strip())
    except ValueError:
        return 0


@click.command()
@click.option('--blastout', help='blast result')
@click.option('--genomefile', help='contigID\tcontiglen')
@click.option('--outfile', help='outfile')
def main(blastout, genomefile, outfile):
    """
    qid col 0
    qstart col 6
    qend col 7
    query id 中不能有括号
    """
    querylen = {x.split()[0]: int(x.split()[1].strip()) for x in open(genomefile).readlines()}
    result = calquerycov(blastout, querylen)
    with open(outfile, 'w') as f:
        for k,v in result.items():
            f.write(f'{k}\t{v}\n')


if __name__ == '__main__':
    main()


