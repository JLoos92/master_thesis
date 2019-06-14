#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 11:14:45 2019

@author: jloos
"""

from __main__ import ModelRun
import matplotlib.pylab as plt
import numpy as np


def compute_correlation_2d(t=None,
                           width=None,
                        **kwargs):

    '''
    Computes root mean square of the hydrostatic deviation. The hydrostatic de-
    viation is saved in an 1d array for the 2d simulation case.       
    '''    
    
    # Import some properties from class    
    list_widths = ModelRun(150,20000,0,'extent',50,"2").list_widths
    list_widths.sort(key=int)
    list_widths = list_widths[0:10]
    original_width = round(np.sqrt(width)*2)*2
    list_widths = [width]
    num_timesteps = ModelRun(150,20000,0,'extent',50,"2").num_timesteps
    num_timesteps = num_timesteps -2
    timesteps = np.arange(0,num_timesteps,1)
    
    
     
    
    # Define rc_params for figure
    fig = plt.figure(figsize = (15,15))        
    params = {
                    'legend.fontsize': 15,
                    'xtick.labelsize': 10,
                    'ytick.labelsize': 10,
                    'axes.labelsize':15,
                    'text.usetex': False,
                    'font.family': 'serif',
                    'axes.titlesize': 18,
                    'axes.titleweight': 'bold'
                    }



    
    for width in list_widths:
        rms_total = []
        rms_total_extent = []
        time = []
        for i in range(2,200):
            
            # Set up class
            mr = ModelRun(150,width,0,0,i,"2")
            mr_extent = ModelRun(150,width,0,'extent',i,'2')
            
            # Calculate hydrostatic thickness
            ht = mr.compute_hydrostatic_thickness()
            ht_extent = mr_extent.compute_hydrostatic_thickness()
            
            
            # compute deviation
            lower = ht[3]
            calc_thickness_bs = ht[1]
            
            lower_extent = ht_extent[3]
            calc_thickness_bs_extent = ht_extent[1]
            
            hydrostatic_deviation = calc_thickness_bs - lower
            hydrostatic_deviation_extent = calc_thickness_bs_extent - lower_extent
            
            rms = np.sqrt(np.mean(hydrostatic_deviation**2))
            rms_extent = np.sqrt(np.mean(hydrostatic_deviation_extent**2))
            
            # Append rms value for timestep i
            rms_total.append(rms)
            rms_total_extent.append(rms_extent)
            
            time.append(i*5)
            
          
            
    #  Make plot
                       
    plt.axes(frameon = 0)
    plt.plot(time,rms_total)
    plt.plot(time,rms_total_extent)
    legend = plt.legend(['cw = ' + str(original_width)+' m', 'cw = ' + str(original_width)+' m (extended domain)'],loc=1)
            
    
            
    frame = legend.get_frame()
    frame.set_facecolor('0.8')
    frame.set_edgecolor('0.8')
    
    
    # Plt properties
    plt.rcParams.update(params) 
    plt.title('RMS of hd with multiple channel widths', y = 1.05)
    plt.xlim(0,1000)
    plt.ylim(0,15)       
    plt.xlabel('Time [a]',labelpad=20)
    plt.ylabel('RMS of hydrostatic deviation [m]',labelpad=20)
    plt.grid(linestyle = '--')
    
    
    path = str('plots/')
    fname= str('corr_2d_' + str(original_width) + '.png')
    
    fig.savefig(path + fname,dpi=300)
    
    
    return rms_total
            
  
    










def compute_maxpeak_dev(t=None,
                        width=None,
                        **kwargs):   
    
    '''
    
    
    
    
    '''
   
    # Define rc_params for figure
    fig = plt.figure(figsize = (15,15))        
    params = {
                    'legend.fontsize': 15,
                    'xtick.labelsize': 10,
                    'ytick.labelsize': 10,
                    'axes.labelsize':15,
                    'text.usetex': False,
                    'font.family': 'serif',
                    'axes.titlesize': 18,
                    'axes.titleweight': 'bold'
                    }

    hd_list_tight = []
    hd_list_wide = []
    time = []

    for i in range(3,200):
        hd_regular = ModelRun(150,width,0,0,i,"2").compute_hydrostatic_thickness()
        hd_extended = ModelRun(150,width,0,'extent',i,"2").compute_hydrostatic_thickness()
        
        ht_regular = hd_regular[7]
        ht_extended = hd_extended[7]
       
        
        hd_list_tight.append(ht_regular)
        hd_list_wide.append(ht_extended)
                
        time.append(i*5)
            
          
            
    #  Make plot
    original_width = round(np.sqrt(width)*2)*2                   
    plt.axes(frameon = 0)
    plt.plot(time,hd_list_tight)
    plt.plot(time,hd_list_wide)
    legend = plt.legend(['cw = ' + str(original_width)+' m', 'cw = ' + str(original_width)+' m (extended domain)'],loc=1)
            
    
            
    frame = legend.get_frame()
    frame.set_facecolor('0.8')
    frame.set_edgecolor('0.8')
    
    
    # Plt properties
    plt.rcParams.update(params) 
    plt.title('The maximum peak deviation', y = 1.05)
    plt.xlim(0,1000)
    plt.ylim(0,100)       
    plt.xlabel('Time [a]',labelpad=20)
    plt.ylabel('Max peak deviation percent',labelpad=20)
    plt.grid(linestyle = '--')
    
    
    path = str('plots/')
    fname= str('maxpeak_dev_2d_' + str(original_width) + '.png')
    
    fig.savefig(path + fname,dpi=300)