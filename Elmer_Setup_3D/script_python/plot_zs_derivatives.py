#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:49:11 2019

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

def plot_zs_derivative(t1 = None,
                     prop = None,
                     x_prof= None):
                
    '''
    Function: plot_zs_derivative
    ---------------
    Additional plot function. Shows in the upper subplot the first derivative
    of the top surface zs. The lower subplot is a line plot of a bridging 
    profile at x_prof position (default is -124 m; x_prof = 7).
    
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
        t1 = 50
        
    # Default is zero for regular domain size
    if prop is None:
        prop=0
    
    # Transect of profile
    if x_prof is None:
        x_prof = 7
         
    
    


    # Colors
    orange = '#D55E00'  
    
    nrow = 2; ncol = 1;
    fig, axsi = plt.subplots(nrow,ncol,sharex = True, squeeze = False, constrained_layout=True)
    
    list_widths = list(np.arange(1000,10000,1000))
    list_widths_1 = list(np.arange(10000,100000,10000))
    list_widths_2 = list(np.arange(100000,900000,100000))
    
    list_widths.extend(list_widths_1)
    list_widths.extend(list_widths_2)
    

    colours = ['black',orange]        
   
    
    for width in list_widths:
        
        original_halfwidth = int(round(np.sqrt(width)*2))
        original_width = original_halfwidth
             
        ax1 = axsi[0][0]
        ax2 = axsi[1][0]
        
        plt.subplots_adjust(hspace=0.1)
        
        
        
        
        # Custom params load from __plot_params
        plt.rcParams.update(params_horizontal) 
        
    
        mr = ModelRun(150,width,0,prop,t1,"2")
        ht = mr.compute_hydrostatic_thickness()
            
        # compute deviation regular
        upper=ht[2]
        points=ht[4]
        x_line = ht[0]
        lower = ht[3]
        calc_thickness_bs = ht[1] 
 
        scalar_real = mr.get_scalar_real(' sxy')
    
        x_mat = scalar_real[0]
        y_mat = scalar_real[1]
        scalar_mat = scalar_real[2]

        
        # Calculated hydrostatic thickness
        calc_thickness_bs = ht[1]
        hydrostatic_deviation = calc_thickness_bs - lower  
        points = ht[4]
        x = points[:,0]
        y = points[:,1]
        
        # place text box in upper left in axes coords      
        props = dict(boxstyle='round', facecolor='wheat')
        ax1.text(0.03, 0.92, str('(a)') + '  t = ' + str(t1*5), transform=ax1.transAxes, 
        verticalalignment='top', bbox=props,weight='bold')
        

        # compute derivatives       
        dx = x_line[1]-x_line[0]
        gf = ndimage.gaussian_filter1d(upper,sigma = 1, order=1, mode='wrap') /dx
        gf_2 = ndimage.gaussian_filter1d(gf,sigma = 1, order=1, mode='wrap') /dx
 
        
        
        ax1.plot(x_line,gf,color = 'black',linestyle = '-',linewidth=0.3)
        ax2.plot(x_mat[x_prof,:],scalar_mat[x_prof,:],linestyle = '-',linewidth = 0.3, color = 'red')
 

        
        
        legend_ax1 = ax1.legend(['Second derivative of $z_{s}$'],loc="upper right", prop=dict(weight='bold'))
        frame_ax1 = legend_ax1.get_frame()
        frame_ax1.set_facecolor('0.7')
        frame_ax1.set_edgecolor('0.7')
        
        legend_ax2 = ax2.legend(['Bridging profile @ -124m'],loc="upper right", prop=dict(weight='bold'))
        frame_ax2 = legend_ax2.get_frame()
        frame_ax2.set_facecolor('0.7')
        frame_ax2.set_edgecolor('0.7')
        
        
        ax1.set_xlim(-2000,2000)
       

        
        ax1.set_ylabel('$\sigma_{xy}$ [MPa]', visible = True)
        ax2.set_ylabel('$\sigma_{xy}$ [MPa]', visible = True)
        ax2.set_xlabel('Across flow distance [m]', visible = True)
        y1,y2 = ax1.get_ylim()
        y_1,y_2 = ax2.get_ylim()
        
               
        fig.add_subplot(ax1)
        fig.add_subplot(ax2)
    
    path = str('plots/06_appendix/')
       
    fname_png = str('zs_firstderivatives_2d_all' + str(t1*5) +'a'+ '.png')
    fname_pdf = str('zs_firstderivatives_2d_all' + str(t1*5) +'a'+'.pdf')

    
    plt.savefig(path + fname_png, format = 'png',dpi=1000,bbox_inches='tight')
    plt.savefig(path + fname_pdf, format = 'pdf',dpi=1000,bbox_inches='tight')       
        
    plt.show()  