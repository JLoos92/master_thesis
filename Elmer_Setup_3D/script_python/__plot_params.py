#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 10:22:09 2019

@author: jloos
"""

import matplotlib.pyplot as plt
import subprocess

subprocess.call(["rsync", "-ruvt","plots","../../../latex_thesis/figures/"])
print('Plot folder synced.')

# plot params
# load plot properties as a module





params = {
                        'legend.fontsize' : 15,
                        'xtick.labelsize': 12,
                        'ytick.labelsize': 12,
                        'axes.labelsize':17,
                        'text.usetex': False,
                        'font.family': 'sans-serif',
                        'axes.titlesize': 20,
                        'axes.titleweight': 'bold',
                        'figure.figsize': [15,15],
                        'figure.frameon':0,
                        'font.sans-serif':'Helvetica',
                        'axes.labelpad':20,
                        'figure.titleweight':'bold',
                        'grid.linestyle': '--',
                        'axes.grid': True
                        
                        }


