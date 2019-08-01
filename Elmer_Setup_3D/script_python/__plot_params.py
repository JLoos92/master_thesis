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
                        'legend.fontsize' : 8,
                        'xtick.labelsize': 9,
                        'xtick.bottom': True,
                        'xtick.top': True,
                        'xtick.labelbottom':True,
                        'xtick.direction':'in',
                        'ytick.direction':'in',
                        'ytick.labelsize': 9,
                        'axes.labelsize':8,
                        'text.usetex': False,
                        'font.family': 'sans-serif',
                        'axes.titlesize': 14,
                        'axes.titleweight': 'bold',
                        'figure.figsize': [4,4],
                        'figure.frameon':0,
                        'font.sans-serif':'Helvetica',
                        'axes.labelpad':10,
                        'figure.titleweight':'bold',
                        'grid.linestyle': '--',
                        'axes.grid': False,
                        'lines.linewidth': 1.5,
                        'xtick.minor.visible': True,
                        'ytick.minor.visible' : True,
                        'figure.autolayout' : False     # False for peak deviation
                        
                        }


params_bridging_2d = {
                        'legend.fontsize' : 4,
                        'xtick.labelsize': 5,
                        'xtick.bottom': True,
                        'xtick.top': True,
                        'xtick.labelbottom':True,
                        'xtick.direction':'in',
                        'ytick.direction':'in',
                        'ytick.labelsize': 5,
                        'axes.labelsize':8,
                        'text.usetex': False,
                        'font.family': 'sans-serif',
                        'axes.titlesize': 14,
                        'axes.titleweight': 'bold',
                        'figure.figsize': [4,4],
                        'figure.frameon':0,
                        'font.sans-serif':'Helvetica',
                        'axes.labelpad':2,
                        'figure.titleweight':'bold',
                        'grid.linestyle': '--',
                        'axes.grid': False,
                        'lines.linewidth': 1,
                        'xtick.minor.visible': False,
                        'ytick.minor.visible' : False,
                        'figure.autolayout' : False     # False for peak deviation
                        
                        }