# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 17:37:56 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import time
#import schedule
import itchat
import os
from shutil import rmtree


def load_dir(dirfile):
    """
    load directory from file
    line started with "#" will be ignored
    """
    dirs = []
    with open(dirfile) as f:
        for line in f:
            if line[0] != '#':
                dirs.append(line.strip())
    return dirs


def clean_dirs(dirs):
    sizes = []
    for directory in dirs:
        print(directory)
        getsize = f"du -sh {directory}"
        print(getsize)
        size = os.popen(getsize).read().split()[0]
        print(size)
        sizes.append(size)
        for file in os.listdir(directory):
            absfile = os.path.join(directory, file)
            if os.path.isdir(absfile):
                rmtree(absfile)
            else:
                os.remove(absfile)
    sent_wechat(dirs, sizes)



def sent_wechat(dirs, sizes):
    content = []
    header = '缓存文件已清理: \n'
    for directory, size in zip(dirs, sizes):
        content.append(f"{size}\t{directory}")
    print(content)
    text = header + '\n'.join(content)
    print(text)
    itchat.send(msg=text, toUserName=itchat.search_friends(name='付玮玮')[0]['UserName'])





@click.command()
@click.option('--dirfile', help='file contains directorys to clean, one dir per line. line started with "#" will be ignored.')
def main(dirfile):
    itchat.auto_login(enableCmdQR=2, hotReload=True)
#    schedule.every().day.at("10:30").do(clean_dirs)
#    schedule.every(0.05).minutes.do(clean_dirs(dirs))
#    while True:
#        schedule.run_pending()
#        time.sleep(1)
    while True:
        dirs = load_dir(dirfile)
        clean_dirs(dirs)
        time.sleep(86400)


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        itchat.logout()
        print(err)