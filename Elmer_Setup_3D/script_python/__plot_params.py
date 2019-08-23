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

# Color palette for plots


params = {
                        'legend.fontsize' : 8,
                        'xtick.labelsize': 9,
                        'xtick.bottom': True,
                        'xtick.top': True,
                        'xtick.labelbottom':True,
                        'xtick.direction':'in',
                        'ytick.direction':'in',
                        'ytick.labelsize': 9,
                        'axes.labelsize':10,
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
                        'legend.fontsize' : 5.7,
                        'xtick.labelsize': 7,
                        'xtick.bottom': True,
                        'xtick.top': False,
                        'xtick.labelbottom':True,
                        'xtick.direction':'in',
                        'ytick.direction':'in',
                        'ytick.labelsize': 7,
                        'axes.labelsize':10,
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
                        'figure.constrained_layout.use' : True     # False for peak deviation
                        
                        }


params_single_2d = {
                        'legend.fontsize' : 9,
                        'xtick.labelsize': 11,
                        'xtick.bottom': True,
                        'xtick.top': False,
                        'xtick.labelbottom':True,
                        'xtick.direction':'in',
                        'ytick.direction':'in',
                        'ytick.labelsize': 11,
                        'axes.labelsize':11,
                        'text.usetex': False,
                        'font.family': 'sans-serif',
                        'axes.titlesize': 14,
                        'axes.titleweight': 'bold',
                        'figure.figsize': [3.7,3.7],
                        'figure.frameon':0,
                        'font.sans-serif':'Helvetica',
                        'axes.labelpad':2,
                        'figure.titleweight':'bold',
                        'grid.linestyle': '--',
                        'axes.grid': False,
                        'lines.linewidth': 1,
                        'xtick.minor.visible': False,
                        'ytick.minor.visible' : False,
                        'figure.constrained_layout.use' : True     # False for peak deviation
                        
                        }



params_3d = {
                        'legend.fontsize' : 6.5,
                        'font.size':5,
                        'font.weight': 'bold',
                        'xtick.labelsize': 5,
                        'xtick.bottom': True,
                        'xtick.top': False,
                        'xtick.labelbottom':True,
                        'xtick.direction':'in',
                        'ytick.direction':'in',
                        'ytick.labelsize': 5,
                        'axes.labelsize':6.5,
                        'text.usetex': True,
                        'axes.titlesize': 14,
                        'axes.titleweight': 'bold',
                        'figure.figsize': [3.7,3.7],
                        'figure.frameon':0,
                        'font.sans-serif':'Helvetica',
                        'axes.labelpad':0.2,
                        'figure.titleweight':'bold',
                        'grid.linestyle': '--',
                        'axes.grid': False,
                        'lines.linewidth': 1,
                        'xtick.minor.visible': False,
                        'ytick.minor.visible' : False,
                        'figure.constrained_layout.use' : True     # False for peak deviation
                        
                        }

params_profile = {
                        'legend.fontsize' : 8,
                        'font.size':5,
                        'font.weight': 'bold',
                        'xtick.labelsize': 7,
                        'xtick.bottom': True,
                        'xtick.top': False,
                        'xtick.labelbottom':True,
                        'xtick.direction':'in',
                        'ytick.direction':'in',
                        'ytick.labelsize': 7,
                        'axes.labelsize':9,
                        'text.usetex': False,
                        'font.family': 'sans-serif',
                        'axes.titlesize': 14,
                        'axes.titleweight': 'bold',
                        'figure.figsize': [6,4],
                        'figure.frameon':0,
                        'font.sans-serif':'Helvetica',
                        'axes.labelpad':1,
                        'figure.titleweight':'bold',
                        'grid.linestyle': '--',
                        'axes.grid': False,
                        'lines.linewidth': 1.5,
                        'xtick.minor.visible': False,
                        'ytick.minor.visible' : False,
                        'figure.constrained_layout.use' : True     # False for peak deviation
                        
                        }

params_horizontal = {
                        'legend.fontsize' : 6,
                        'font.size':6.5,
                        'font.weight': 'bold',
                        'xtick.labelsize': 6.5,
                        'xtick.bottom': True,
                        'xtick.top': False,
                        'xtick.labelbottom':True,
                        'xtick.direction':'in',
                        'ytick.direction':'in',
                        'ytick.labelsize': 6.5,
                        'axes.labelsize':6.5,
                        'text.usetex': True,
                        'axes.titlesize': 6.5,
                        'axes.titleweight': 'bold',
                        'figure.figsize': [3.7,2.4],
                        'figure.frameon':0,
                        'font.sans-serif':'Helvetica',
                        'axes.labelpad':1,
                        'figure.titleweight':'bold',
                        'figure.titlesize':6.5,
                        'grid.linestyle': '--',
                        'axes.grid': False,
                        'lines.linewidth': 1,
                        'xtick.minor.visible': False,
                        'ytick.minor.visible' : False,
                        'figure.constrained_layout.use' : True                      
                        }

params_horizontal_single = {
                        'legend.fontsize' : 9,
                        'font.size':11,
                        'font.weight': 'bold',
                        'xtick.labelsize': 11,
                        'xtick.bottom': True,
                        'xtick.top': False,
                        'xtick.labelbottom':True,
                        'xtick.direction':'in',
                        'ytick.direction':'in',
                        'ytick.labelsize': 11,
                        'axes.labelsize':11,
                        'text.usetex': True,
                        'font.sans-serif': 'Helvetica',
                        'axes.titlesize': 12,
                        'axes.titleweight': 'bold',
                        'figure.figsize': [3.7,3.7],
                        'figure.frameon':0,
                        'font.sans-serif':'Helvetica',
                        'axes.labelpad':1,
                        'figure.titleweight':'bold',
                        'figure.titlesize':11,
                        'grid.linestyle': '--',
                        'axes.grid': False,
                        'lines.linewidth': 2,
                        'xtick.minor.visible': False,
                        'ytick.minor.visible' : False,
                        'figure.constrained_layout.use' : True
                         # False for peak deviation                        
                        }

params_vertical = {
                        'legend.fontsize' : 6.5,
                        'font.size':11,
                        'font.weight': 'bold',
                        'xtick.labelsize': 11,
                        'xtick.bottom': True,
                        'xtick.top': False,
                        'xtick.labelbottom':True,
                        'xtick.direction':'in',
                        'ytick.direction':'in',
                        'ytick.labelsize': 11,
                        'axes.labelsize':11,
                        'text.usetex': True,
                        'font.sans-serif': 'Helvetica',
                        'axes.titlesize': 12,
                        'axes.titleweight': 'bold',
                        'figure.figsize': [6,3],
                        'figure.frameon':0,
                        'font.sans-serif':'Helvetica',
                        'axes.labelpad':1,
                        'figure.titleweight':'bold',
                        'figure.titlesize':11,
                        'grid.linestyle': '--',
                        'axes.grid': False,
                        'lines.linewidth': 1.5,
                        'xtick.minor.visible': False,
                        'ytick.minor.visible' : False,
                        'figure.constrained_layout.use' : True
                         # False for peak deviation                        
                        }