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
                        'legend.fontsize' : 12,
                        'xtick.labelsize': 20,
                        'xtick.bottom': True,
                        'xtick.labelbottom':True,
                        'xtick.direction':'in',
                        'ytick.labelsize': 20,
                        'axes.labelsize':25,
                        'text.usetex': False,
                        'font.family': 'sans-serif',
                        'axes.titlesize': 30,
                        'axes.titleweight': 'bold',
                        'figure.figsize': [100,100],
                        'figure.frameon':0,
                        'font.sans-serif':'Helvetica',
                        'axes.labelpad':18,
                        'figure.titleweight':'bold',
                        'grid.linestyle': '--',
                        'axes.grid': True,
                        'lines.linewidth': 6
                        
                        }


