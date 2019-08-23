#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 14:02:38 2019

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
import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

    
def profile_velos_2d(t1 = None,
                     width_1 = None,
                     width_2 = None,
                     prop = None,
                     x_prof= None):
                

    
    
    
        # Choose default timesteps if timesteps are not given
    if t1 is None:
        t1 = 50

     
    if prop is None:
        prop=0
    
    if x_prof is None:
        x_prof = 7
         
    



    # Colors
    orange = '#D55E00'  
    red = '#CC2529'
    wheat = '#948B3D'
    
    nrow = 2; ncol = 1;
    fig, axsi = plt.subplots(nrow,ncol,sharex = True, squeeze = False, constrained_layout=True)
    
    
    widths = [width_1,width_2]  
    colours = ['black',orange]        
   
    
    for width,colour in zip(widths,colours):
        
        original_halfwidth = int(round(np.sqrt(width)*2))
        original_width = original_halfwidth
        
        original_halfwidth_1 = int(round(np.sqrt(width_1)*2))
        original_halfwidth_2 = int(round(np.sqrt(width_2)*2))
    
       
        
        ax1 = axsi[0][0]
        ax2 = axsi[1][0]
        
        plt.subplots_adjust(hspace=0.1)
        
        
        
        
         # Custom params load from __plot_params
        plt.rcParams.update(params_horizontal) 
        
    
        mr = ModelRun(150,width,0,prop,t1,"2")
        ht = mr.compute_hydrostatic_thickness()
    
        upper=ht[2]
        lower=ht[3]
        
        scalar_real = mr.get_scalar_real('velocity')
        
        x_mat = scalar_real[0]
        y_mat = scalar_real[1]
        scalar_mat_x = scalar_real[2]
        scalar_mat_y = scalar_real[3]
        
        
        # Calculated hydrostatic thickness
        calc_thickness_bs = ht[1]
          
        points = ht[4]
        x = points[:,0]
        y = points[:,1]
        
        # place text box in upper left in axes coords      
        props = dict(boxstyle='round', facecolor='wheat')
        ax1.text(0.03, 0.92, str('(a)') + '  t = ' + str(t1*5) + 'a\n @ y =  ' + str(int(np.mean(y_mat[x_prof,:]))) + 'm', transform=ax1.transAxes, 
        verticalalignment='top', bbox=props,weight='bold')
        
        ax2.text(0.03, 0.92, str('(b)') + '  t = ' + str(t1*5) + 'a\n @ y =  ' + str(int(np.mean(y_mat[x_prof,:]))) + 'm', transform=ax2.transAxes,
        verticalalignment='top', bbox=props,weight='bold')
       
        
        ax1.plot(x_mat[x_prof,:],scalar_mat_x[x_prof,:],color = colour,linestyle = '-',linewidth=1)
        ax2.plot(x_mat[x_prof,:],scalar_mat_y[x_prof,:],color = colour,linestyle = '-',linewidth=1)
        
        legend_ax1 = ax1.legend(['cw = ' + str(original_halfwidth_1*2) + 'm','cw = ' + str(original_halfwidth_2*2) + 'm'],loc="upper right", prop=dict(weight='bold'))
        frame_ax1 = legend_ax1.get_frame()
        frame_ax1.set_facecolor('0.7')
        frame_ax1.set_edgecolor('0.7')
        
        
        ax1.set_xlim(-2500,2500)
        
        ax1.set_ylabel('$u_{x}$ [$m\: s^{-1}$]', visible = True)
        ax2.set_ylabel('$u_{y}$ [$m\: s^{-1}$]', visible = True)
        ax2.set_xlabel('Across flow distance [m]', visible = True)
        y1,y2 = ax1.get_ylim()
        y_1,y_2 = ax2.get_ylim()
        
    
            
        ax1.vlines(-original_width,y1,y2,linewidth=0.5,color = red)
        ax1.vlines(original_width,y1,y2,linewidth=0.5,color = red)
        ax2.vlines(-original_width,y_1,y_2,linewidth=0.5,color = red)
        ax2.vlines(original_width,y_1,y_2,linewidth=0.5,color = red)

   
    # Save figures    
    path = str('plots/Final_plots/')
       
    fname_png = str('profile_velos_2d_' + str(original_width*2) + '.png')
    fname_pdf = str('profile_velos_2d_' + str(original_width*2) + '.pdf')

    
    plt.savefig(path + fname_png, format = 'png',dpi=1000,bbox_inches='tight')
    plt.savefig(path + fname_pdf, format = 'pdf',dpi=1000,bbox_inches='tight')       
        
    plt.show()    