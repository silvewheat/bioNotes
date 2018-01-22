# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 22:13:34 2017

@author: Caiyd
"""

import os
import click
import gzip
import pandas as pd
from itertools import combinations
from collections import OrderedDict




@click.command()
@click.argument('listfile')
@click.argument('outfile')
@click.option('--n-col', help='column to use, 0 based', type=int, default=0)
@click.option('--outsites', help='file prefix to save shared sites, default is False', default=False)
def main(listfile, outfile, n_col, outsites):
    """
    listfile 是文件列表，一行一个文件，如果不在当前目录的话需要加上路径。里面的文件默认是按空白符分列，文件名不能有重复
    outfile 是输出的统计文件
    """
    filelist = [x.strip() for x in open(listfile).readlines()]
    set_dict = OrderedDict()
    idlist = []
    for file in filelist:
        id_ = os.path.basename(file)
        idlist.append(id_)
        set_dict[id_] = set(pd.read_csv(file, sep='\t', usecols=[n_col], squeeze=True, dtype=str, header=None).tolist())
    with open(outfile, 'w') as f:
        for n_pair in range(1, len(filelist)+1):
            for pair in combinations(idlist, n_pair):
                tmp_set = set()
                tmp_set |= set_dict[pair[0]]
                for item in pair:
                    tmp_set &= set_dict[item]
                n_shared = len(tmp_set)
                f.write('_'.join(pair) + '\t%s\n' % n_shared)
                if outsites and (len(pair)>1):
                    with gzip.open(f'{outsites}_' + '_'.join(pair) + '.sites.gz', 'wb') as f2:
                        f2.write('\n'.join(tmp_set).encode())


if __name__ == '__main__':
    main()
