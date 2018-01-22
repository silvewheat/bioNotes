# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 22:38:13 2018

@author: Caiyd
"""

import os
import sys
import time
import click
import numpy as np



def check_jobs(jobids):
    """
    返回没有DONE的任务数量
    返回值为0意味着都DONE了
    要是有EXIT的就重新提交那个任务
    """
    stats = []
    exists = []
    newjobs = []
    for njob, jobid in enumerate(jobids, 1):
        cmd = 'jjobs %s' % jobid
        stat = os.popen(cmd).readlines()[1].split()[2]
        if stat != 'EXIT':
            stats.append(stat)
        else:
            exists.append(jobid)
            cmd = f'jsub -R "rusage[res=1]span[hosts=1]" -q jynodequeue -e 01.shell/dadi.e -o 01.shell/dadi.o -n 2 -J dadi_{njob} -M 20000000 bash 01.shell/{njob}.sh'
            newid = os.popen(cmd).readlines()[-1].strip().split()[1][1:-1]
            newjobs.append(newid)
            stats.append('RUN')
            print(f'{jobid} was EXIST, and resubmitted({newid})')
    if exists:
        for jobid in exists:
            jobids.remove(jobid)
        jobids.extend(newjobs)
    return np.sum(np.array(stats) != 'DONE'), jobids


def produce_shell(i):
    return f"""source /stor9000/apps/users/NWSUAF/2015060152/.bashrc
/stor9000/apps/appsoftware/BioSoftware/bin/python2 \\
    {os.getcwd()}/01.model.py \\
    {os.getcwd()}/dadi.fs \\
    > {os.getcwd()}/output_dadi/output.{i}.txt """



def submit():
    jobids = []
    for i in range(1, 101):
        cmd = f'jsub -R "rusage[res=1]span[hosts=1]" -q jynodequeue -e 01.shell/dadi.e -o 01.shell/dadi.o -n 2 -J dadi_{i} -M 20000000 bash 01.shell/{i}.sh'
        line = os.popen(cmd).readlines()[-1]
        jobids.append(line.strip().split()[1][1:-1])
    return jobids


@click.command()
def main():
    """
    有程序exit就死循环了
    """
    os.mkdir(r'output_dadi')
    os.mkdir(r'01.shell')
    for i in range(1, 101):
        with open(f'01.shell/{i}.sh', 'w') as f:
            f.write(produce_shell(i))
    for i in range(1, 4):
        print(f'batch{i}')
        jobids = submit()
        print('submitted:')
        print(' '.join(jobids))
        running = True
        while running:
            n_running, jobids = check_jobs(jobids)
            sys.stdout.write(' ' * 30 + '\r')
            sys.stdout.flush()
            sys.stdout.write(f'{n_running} jobs running' + '\r')
            sys.stdout.flush()
            if n_running == 0:
                running = False
            else:
                time.sleep(1)
        os.system('perl 03.collect.pl')

if __name__ == '__main__':
    main()




