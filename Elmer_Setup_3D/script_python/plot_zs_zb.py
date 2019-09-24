#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 21:06:22 2019

@author: jloos
"""
from main import ModelRun
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from scipy.interpolate import griddata
import pandas as pd
from __plot_params import params_horizontal
from scipy.interpolate import interp1d
from numpy import diff
from scipy import ndimage
import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

    
def plot_zs_zb(t1 = None,
              prop = None):
                
                
    '''
    Function: plot_zs_zb
    ---------------
    Two subplots. Upper subplots contains all top surfaces zs. The lower sub-
    plot shows all lower surfaces zb. Timestep can be chosen.
    
    Parameters
    ----------
    t1 : int
        timestep for plot
    prop : int (0 for regular domain)
        timestep for plot  
    x_prof : int
        transect of the ice-shelf for the 'n'th column (15 layers)
        default is 7, which is approximately -124 m, closely above the channel 
        apex
        
    '''
    
    
    
    # Choose default timesteps if timesteps are not given
    if t1 is None:
        t1 = 40
     
    if prop is None:
        prop=0
    
    
    # Colors
    orange = '#D55E00'  
    red = '#CC2529'
    wheat = '#948B3D'
    
    nrow = 2; ncol = 1;
    fig, axsi = plt.subplots(nrow,ncol,sharex = True, squeeze = False, constrained_layout=True)
    
    list_widths = list(np.arange(1000,10000,1000))
    list_widths_1 = list(np.arange(10000,100000,10000))
    list_widths_2 = list(np.arange(100000,900000,100000))
    
    list_widths.extend(list_widths_1)
    list_widths.extend(list_widths_2)

    for width in list_widths:
        
        original_halfwidth = int(round(np.sqrt(width)*2))
 
        # define subplots
        ax1 = axsi[0][0]
        ax2 = axsi[1][0]
        
        # vertical space between subplots
        plt.subplots_adjust(hspace=0.1)
       
        # Custom params load from __plot_params
        plt.rcParams.update(params_horizontal) 
        
    
        mr = ModelRun(150,width,0,prop,t1,"2")
        ht = mr.compute_hydrostatic_thickness()
            
        # compute deviation regular
        upper=ht[2]
        x_line = ht[0]
        lower = ht[3]

        
        # place text box in upper left in axes coords      
        props = dict(boxstyle='round', facecolor='wheat')
        ax1.text(0.03, 0.92, str('(b)') + '  t = ' + str(t1*5), transform=ax1.transAxes, 
        verticalalignment='top', bbox=props,weight='bold')


        # Plot upper and lower surface
        ax1.plot(x_line,upper,color = 'black',linestyle = '-',linewidth=0.3)
        ax2.plot(x_line,lower,linestyle = '-',linewidth = 0.3, color = 'black')
 

        # compute derivatives       
        dx = x_line[1]-x_line[0]
        gf = ndimage.gaussian_filter1d(upper,sigma = 5, order=1, mode='wrap') /dx
        gf_2 = ndimage.gaussian_filter1d(gf,sigma = 1, order=1, mode='wrap') /dx




        # upper and lower legend
        legend_ax1 = ax1.legend(['$z_{s}$'],loc="upper right", prop=dict(weight='bold'))
        frame_ax1 = legend_ax1.get_frame()
        frame_ax1.set_facecolor('0.7')
        frame_ax1.set_edgecolor('0.7')
        
        legend_ax2 = ax2.legend(['$z_{b}$'],loc="upper right", prop=dict(weight='bold'))
        frame_ax2 = legend_ax2.get_frame()
        frame_ax2.set_facecolor('0.7')
        frame_ax2.set_edgecolor('0.7')
        
        
        # set lims
        ax1.set_xlim(-5000,5000)
       

        # label and get lims
        ax1.set_ylabel('Shelf elevation [m]', visible = True)
        ax2.set_ylabel('Shelf elevation [m]', visible = True)
        ax2.set_xlabel('Across flow distance [m]', visible = True)
        y1,y2 = ax1.get_ylim()
        y_1,y_2 = ax2.get_ylim()
        
               
        fig.add_subplot(ax1)
        fig.add_subplot(ax2)
    
    path = str('plots/04_discussion/')
       
    fname_png = str('zs_derivatives_2d_' + str(t1*5) +'a'+ '.png')
    fname_pdf = str('zs_derivatives_2d_' + str(t1*5) +'a'+'.pdf')

    
    plt.savefig(path + fname_png, format = 'png',dpi=1000,bbox_inches='tight')
    plt.savefig(path + fname_pdf, format = 'pdf',dpi=1000,bbox_inches='tight')       
        
    plt.show() 
