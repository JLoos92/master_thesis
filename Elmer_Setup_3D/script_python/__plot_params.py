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
                        'legend.fontsize' : 6,
                        'xtick.labelsize': 14,
                        'xtick.bottom': True,
                        'xtick.top': True,
                        'xtick.labelbottom':True,
                        'xtick.direction':'in',
                        'ytick.direction':'in',
                        'ytick.labelsize': 14,
                        'axes.labelsize':14,
                        'text.usetex': False,
                        'font.family': 'sans-serif',
                        'axes.titlesize': 30,
                        'axes.titleweight': 'bold',
                        'figure.figsize': [6,6],
                        'figure.frameon':0,
                        'font.sans-serif':'Helvetica',
                        'axes.labelpad':18,
                        'figure.titleweight':'bold',
                        'grid.linestyle': '--',
                        'axes.grid': False,
                        'lines.linewidth': 2,
                        'xtick.minor.visible': True,
                        'ytick.minor.visible' : True,
                        'figure.autolayout' : True
                        
                        }


