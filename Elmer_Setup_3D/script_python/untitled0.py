#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:33:32 2019

@author: jloos
"""

from __main__ import ModelRun
import pandas as pd
from scipy.interpolate import griddata

from __main__ import ModelRun
from __plot_params import params

#
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd





def compute_total_hd_2d(t=None,
                        **kwargs):

    '''
    Computes root mean square of the hydrostatic deviation. The hydrostatic de-
    viation is saved in an 1d array for the 2d simulation case.       
    '''    
    
    # Import some properties from class    
    list_widths = ModelRun(150,20000,0,'extent',50,"2").list_widths
    list_widths.sort(key=int)
    list_widths = list_widths[0:10]
    
    # List of widths for extended and tight domain
    list_widths = list(np.arange(5000,10000,1000))
    list_widths_1 = list(np.arange(20000,100000,10000))
    list_widths_2 = list(np.arange(200000,1000000,100000))
    
    list_widths.extend(list_widths_1)
    
    list_widths.extend(list_widths_2)
    
    print(list_widths)
    
    num_timesteps = ModelRun(150,20000,0,'extent',50,"2").num_timesteps
    num_timesteps = num_timesteps -2
   
    
    times = [5,25,50,100,150]
    
    # Define rc_params for figure
    fig, ax1 = plt.subplots()  
    
    for i in times:
        
        
        
        rms_total = []
        rms_total_extent = []
        time = []
        widths_new = []
    
    
        for width in list_widths:
            
            # Set up class
            mr = ModelRun(150,width,0,0,i,"2")
            #mr_extent = ModelRun(150,width,0,'extent',i,'2')
            
            # Calculate hydrostatic thickness
            ht = mr.compute_hydrostatic_thickness()
            #ht_extent = mr_extent.compute_hydrostatic_thickness()
            
            
            # compute deviation
            lower = ht[3]
            calc_thickness_bs = ht[1]
            
            #lower_extent = ht_extent[3]
            #calc_thickness_bs_extent = ht_extent[1]
            
            hydrostatic_deviation = calc_thickness_bs - lower
            #hydrostatic_deviation_extent = calc_thickness_bs_extent - lower_extent
            
            # Adapt ratio of wider domain
            #hydrostatic_deviation_extent = hydrostatic_deviation_extent[166:501]
            
            rms = np.sqrt(np.mean(hydrostatic_deviation**2))
            #rms_extent = np.sqrt(np.mean(hydrostatic_deviation_extent**2))
            
            # Append rms value for timestep i
            rms_total.append(rms)
            #rms_total_extent.append(rms_extent)
            original_width = int(round(np.sqrt(width)*2)*2)
            widths_new.append(original_width)
            
          
            
        # Make plot
                           
        #ax1.axes(frameon = 0)
        ax1.plot(widths_new,rms_total,linewidth=4)
        legend = ax1.legend(['t = ' + str(5*5)+' a','t = ' + str(25*5)+' a','t = ' + str(50*5)+' a', 't = ' + str(100*5)+' a','t = ' + str(150*5)+' a'],loc=1)
                
        
                
        frame = legend.get_frame()
        frame.set_facecolor('0.7')
        frame.set_edgecolor('0.7')
        
        
        
        
        # Plt properties
        plt.rcParams.update(params) 
        # Plt properties
        plt.rcParams.update(params)
        ax1.spines['top'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(direction='in',length=6,width=2)
        
        #plt.title('RMS of hd with multiple channel widths', y = 1.05)
        ax1.set_ylim(0,8)       
        ax1.set_xlabel('Channel widths' + ' [m]',labelpad=20)
        ax1.set_ylabel('RMS of hydrostatic deviation [m]',labelpad=20)
        ax1.grid(linestyle = '--')
        
        
        path = str('plots/')
        fname= str('allwidths_corr_2d' + '.eps')
        
        fig.savefig(path + fname, format = 'eps', dpi=1000)
    
