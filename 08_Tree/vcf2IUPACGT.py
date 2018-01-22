# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 21:01:14 2018

@author: Caiyd
"""

import os
import click



def bcftools_query():
    return """ bcftools query -f '[ %IUPACGT]\n' """