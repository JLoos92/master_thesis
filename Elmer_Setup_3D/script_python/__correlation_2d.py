#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 11:14:45 2019

@author: jloos
"""

from __main__ import ModelRun
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import pandas as pd

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
    
    # List of widths for extended and tight domain
    list_widths = list(np.arange(10000,100000,10000))
    
    
    
    num_timesteps = ModelRun(150,20000,0,'extent',50,"2").num_timesteps
    num_timesteps = num_timesteps -2
    timesteps = np.arange(0,num_timesteps,1)
    
    
    for width in list_widths:
        rms_total = []
        rms_total_extent = []
        time = []
        original_width = round(np.sqrt(width)*2)*2
     
    
        # Define rc_params for figure
        fig = plt.figure(figsize = (15,15))        
        params = {
                    'legend.fontsize': 15,
                    'xtick.labelsize': 11,
                    'ytick.labelsize': 11,
                    'axes.labelsize':16,
                    
                    'font.family': 'serif',
                    'axes.titlesize': 18,
                    'axes.titleweight': 'bold'
                    }



    
    
    
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
            
            # Adapt ratio of wider domain
            hydrostatic_deviation_extent = hydrostatic_deviation_extent[166:501]
            
            rms = np.sqrt(np.mean(hydrostatic_deviation**2))
            rms_extent = np.sqrt(np.mean(hydrostatic_deviation_extent**2))
            
            # Append rms value for timestep i
            rms_total.append(rms)
            rms_total_extent.append(rms_extent)
            
            time.append(i*5)
            
          
            
        # Make plot
                           
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
   
    subprocess.call(["rsync", "-ruvt","plots/*","../../../latex_thesis/figures/"])
    
    points = ModelRun(150,20000,0,0,50,"2").Points
    points_extent =  ModelRun(150,20000,0,'extent',50,"2").Points
    
    x = points[:,0]
    x_extent = points_extent[:,0]
    
     # List of widths for extended and tight domain
    list_widths = list(np.arange(10000,100000,10000))
    
    for width in list_widths:
        # Define rc_params for figure
        fig, ax1 = plt.subplots()        
        params = {
                        'legend.fontsize': 15,
                        'xtick.labelsize': 12,
                        'ytick.labelsize': 12,
                        'axes.labelsize':17,
                        'text.usetex': False,
                        'font.family': 'serif',
                        'axes.titlesize': 18,
                        'axes.titleweight': 'bold',
                        'figure.figsize': [15,15],
                        'figure.frameon':0
                        
                        }
    
        hd_list_tight = []
        hd_list_wide = []
        sxy_list = []
        sxy_list_extent = []
        time = []

        for i in range(3,200):
            hd_regular = ModelRun(150,width,0,0,i,"2").compute_hydrostatic_thickness()
            hd_extended = ModelRun(150,width,0,'extent',i,"2").compute_hydrostatic_thickness()
            
            # Get bridging stresses
            s_xy =  ModelRun(150,width,0,0,i,"2").sxy
            s_xy_extent = ModelRun(150,width,0,'extent',i,"2").sxy
            
            df = pd.DataFrame({'sxy':s_xy,'x':x})
            df = df[df.x>0]
            sxy_array = df.iloc[:,0].values
            rms_sxy = np.sqrt(np.mean(sxy_array**2))
            
            df_extent = pd.DataFrame({'sxy_extent':s_xy_extent,'x_extent':x_extent})
            df_extent = df_extent[df_extent.x_extent>0]
            df_extent = df_extent[df_extent.x_extent<10000]
            sxy_array_extent = df_extent.iloc[:,0].values
            rms_sxy_extent = np.sqrt(np.mean(sxy_array_extent**2))
            
            
            # Get values of max peak deviation (channel's amplitude) in percent
            ht_regular = hd_regular[7]
            ht_extended = hd_extended[7]
            
            # Append values
            hd_list_tight.append(ht_regular)
            hd_list_wide.append(ht_extended)
            
            sxy_list.append(rms_sxy)
            sxy_list_extent.append(rms_sxy_extent)
                    
            time.append(i*5)
            
          
            
        #  Make plot, cut frame, properties defined in rc.params
        original_width = int(round(np.sqrt(width)*2)*2)                   
        #ax1.set_frame_on(frameon = False)
        ax1.plot(time,hd_list_tight,'b-')
        ax1.plot(time,hd_list_wide,'r-')
        legend = ax1.legend(['cw = ' + str(original_width)+' m (regular domain)', 'cw = ' + str(original_width)+' m (extended domain)'],loc=1)
                
        
                
        frame = legend.get_frame()
        frame.set_facecolor('0.8')
        frame.set_edgecolor('0.8')
        
        
        # Plt properties for ax1 and peak dev
        plt.rcParams.update(params) 
        #ax1.title('The maximum peak deviation', y = 1.05)
        ax1.set_xlim(0,1000)
        ax1.set_ylim(0,100)       
        ax1.set_xlabel('Time [a]',labelpad=20)
        ax1.set_ylabel('Max peak deviation [%]',labelpad=20)
        ax1.grid(linestyle = '--')
        
        # Plt properties for ax2 (bridging)  
        ax2 = ax1.twinx()
        
        ax2.plot(time, sxy_list, 'b--')
        ax2.plot(time, sxy_list_extent, 'r--')
        ax2.set_ylabel('$\sigma_{xy}$ [MPa]',labelpad = 20)
        ax2.tick_params('y')
        
        
        legend_ax2 = ax2.legend(['Bridging (regular domain)','Bridging (extended domain)'],loc=2)
        frame_ax2 = legend_ax2.get_frame()
        frame_ax2.set_facecolor('0.8')
        frame_ax2.set_edgecolor('0.8')
        
        plt.title('Peak deviation and bridging stresses for a channel width of ' + str(original_width) + ' m', y = 1.05)
        
        path = str('plots/')
        fname= str('maxpeak_dev_2d_' + str(original_width) + '.eps')
        
        fig.savefig(path + fname,format = 'eps',dpi=1000)
    
    
    
    
    
    
    
    
    
    