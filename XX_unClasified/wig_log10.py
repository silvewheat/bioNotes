# -*- coding: utf-8 -*-
"""
Created on Tue May 29 22:36:47 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np



@click.command()
@click.argument('wigfile')
def main(wigfile):
    with open(wigfile) as f:
        for line in f:
            if line[0] != '#':
                tline = line.strip().split()
                value = float(tline[-1])
                if value >= 1:
                    tline[-1] = str(np.log10(value))
                    print('\t'.join(tline))
                else:
                    print(line, end='')
            else:
                print(line, end='')

if __name__ == '__main__':
    main()
