#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 17:08:06 2019

@author: jloos
"""

# Custom modules
from main import ModelRun
from __plot_params import params_horizontal
from local_minimum import local_minimum_1d
#
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'





def compute_abutment_peakdeviation(t_end=None,
                                   **kwargs):   
    
    '''
      
    Function: compute_maxpeak_dev 
    
    Parameters
    ----------
    t : int
        timestep        
    
    width : kwarg
        width of a channel
      
    
    '''
   
   
    # Get domain-size
    points = ModelRun(150,20000,0,0,50,"2").Points

    x = points[:,0]
    width_1 = 50000
    width_2 = 100000
    # Time + numbers for ax
    nums = ["(a)","(b)","(c)","(d)","(e)","(f)"]
    
    orange = '#D55E00'
    fig, axs = plt.subplots(1,1, sharex=True, sharey = True) 
    time = [t_end]
    
     # Custom model load from __plot_params
    plt.rcParams.update(params_horizontal)
    plt.subplots_adjust(hspace=0,wspace=0)
    
    abutment_narrow = []
    abutment_wide = []
   
    peak_narrow = []
    peak_wide = []
    
    times = []
    
    #  Make plot, cut frame, properties defined in rc.params
    def original_width(width):
        original_width = int(round(np.sqrt(width)*2)*2)
        return original_width      
    
    
    for t in range(3,t_end):
        
        mr_narrow = ModelRun(150,width_1,0,0,t,"2")
        mr_wide = ModelRun(150,width_2,0,0,t,"2")
        
        ht_narrow =mr_narrow.compute_hydrostatic_thickness()
        ht_wide = mr_wide.compute_hydrostatic_thickness()
        
        hd_n = ht_narrow[7]
        hd_w = ht_wide[7]       
                
        peak_narrow.append(hd_n)
        peak_wide.append(hd_w)
        
        abutment_narrow.append(ht_narrow[9])
        abutment_wide.append(ht_wide[9])
        
        times.append(t*5)
        
    
    
    
    axs.plot(times,peak_narrow,'k-')
    axs.plot(times,peak_wide,'k:')
    axs.set_ylim(0,22)
    
    axs.set_ylabel('Channel peak deviation [\%]')    
    axs.set_xlabel('Time [a]')    
    legend = axs.legend(['cw = ' + str(original_width(width_1)),'cw = ' + str(original_width(width_2))],loc='upper right')                
    frame = legend.get_frame()
    frame.set_facecolor('0.7')
    frame.set_edgecolor('0.7')   
    
    ax1 = axs.twinx()
    ax1.plot(times,abutment_narrow,'r-')
    ax1.plot(times,abutment_wide,'r:')
    ax1.set_ylim(0,22)
    ax1.set_ylabel('Channel abutment deviation [\%]')   
    
        
      
    
    path = str('plots/Final_plots/')
    fname_eps = str('peak_abutment_2d_' + str(width_1)+ str(width_2) + '_' + str(t_end*5) + 'a' + '.pdf')
    fname_png = str('peak_abutment_2d_'+ str(width_1)+ str(width_2) + '_' + str(t_end*5) + 'a' + '.png')
    
    fig.savefig(path + fname_eps, format = 'pdf',dpi=1000,bbox_inches = 'tight') 
    fig.savefig(path + fname_png, format = 'png',dpi=1000,bbox_inches = 'tight') 
    
    
    plt.show()