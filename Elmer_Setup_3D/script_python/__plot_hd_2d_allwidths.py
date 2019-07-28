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





def compute_allwidths_hd_2d(t=None,
                        **kwargs):

    '''
    Computes rms of the hydrostatic deviation. The hydrostatic de-
    viation is saved in an 1d array for the 2d simulation case.       
    '''    
    
    
    # List of widths for extended and tight domain
    list_widths = list(np.arange(5000,10000,1000))
    list_widths_1 = list(np.arange(10000,100000,10000))
    list_widths_2 = list(np.arange(200000,1000000,100000))
    
    list_widths.extend(list_widths_1)
    
    list_widths.extend(list_widths_2)
    
    print(list_widths)
    
    num_timesteps = ModelRun(150,20000,0,'extent',50,"2").num_timesteps
    num_timesteps = num_timesteps -2
   
    
    times = [t]
    
    # Define rc_params for figure
    fig, ax1 = plt.subplots()  
    
    for i in times:
        
        
        # Hydrostatic deviation
        rms_total = []
        rms_total_extent = []
        
        
        # Velocities
        rms_uy_normal = []
        rms_uy_extent = []
        rms_uy_wideextent = []
        
        
        time = []
        widths_new = []
    
    
        for width in list_widths_1:
            
            
            # Define width in m
            original_width = int(round(np.sqrt(width)*2))
            widths_new.append(original_width)
            
            # Set up class
            mr = ModelRun(150,width,0,0,i,"2")
            mr_extent = ModelRun(150,width,0,'extent',i,'2')
            mr_wideextent = ModelRun(150,width,0,'wideextent',i,'2')
            
            # Get velocities for y
            uy_normal = mr.get_scalar('velocity')
            uy_normal_y = uy_normal[:,1] 
            uy_extent = mr_extent.get_scalar('velocity')
            uy_extent_y = uy_extent[:,1] 
            uy_wideextent = mr_wideextent.get_scalar('velocity')
            uy_wideextent_y = uy_wideextent[:,1] 
            
            
            #Get points
            points_normal = mr.Points
            points_extent =  mr_extent.Points
            points_wideextent =  mr_wideextent.Points 
            x_normal = points_normal[:,0]
            x_extent = points_extent[:,0]
            x_wideextent = points_wideextent[:,0]
            
            
            # Normal uy
            df_normal = pd.DataFrame({'x_normal':x_normal,'uy_normal':uy_normal_y})
            df_normal = df_normal[df_normal.x_normal<original_width]
            df_normal = df_normal[df_normal.x_normal>-original_width]
            uy_array_normal = df_normal.iloc[:,1].values
            uy_array_normal = np.sqrt(np.mean(uy_array_normal**2))
            
            # Extended uy
            df_extent = pd.DataFrame({'x_extent':x_extent,'uy_extent':uy_extent_y})
            df_extent = df_extent[df_extent.x_extent<original_width]
            df_extent = df_extent[df_extent.x_extent>-original_width]
            uy_array_extent = df_extent.iloc[:,1].values
            uy_array_extent = np.sqrt(np.mean(uy_array_extent**2))
            
            # Wide extended uy
            df_wideextent = pd.DataFrame({'x_wideextent':x_wideextent,'uy_extent':uy_wideextent_y})
            df_wideextent = df_wideextent[df_wideextent.x_wideextent<original_width]
            df_wideextent = df_wideextent[df_wideextent.x_wideextent>-original_width]
            uy_array_wideextent = df_wideextent.iloc[:,1].values
            uy_array_wideextent = np.sqrt(np.mean(uy_array_wideextent**2))
            
                        
            # Calculate hydrostatic thickness
            ht = mr.compute_hydrostatic_thickness()
            ht_extent = mr_extent.compute_hydrostatic_thickness()
            
            
            
            # compute deviation reg7ular
            lower = ht[3]
            calc_thickness_bs = ht[1]
            
            #extent
            lower_extent = ht_extent[3]
            calc_thickness_bs_extent = ht_extent[1]

            
            hydrostatic_deviation = calc_thickness_bs - lower
            hydrostatic_deviation_extent = calc_thickness_bs_extent - lower_extent
        
            
            # Adapt ratio of wider domain
            hydrostatic_deviation_extent = hydrostatic_deviation_extent[166:501]
            
            
            
            
            # Adapt ratio of extra-wide domain
        
            rms = np.sqrt(np.mean(hydrostatic_deviation**2))
            rms_extent = np.sqrt(np.mean(hydrostatic_deviation_extent**2))
            
           
            # Append rms value for timestep i
            rms_total.append(rms)
            rms_total_extent.append(rms_extent)
            
            # Append uy values
            rms_uy_normal.append(uy_array_normal)
            rms_uy_extent.append(uy_array_extent)
            rms_uy_wideextent.append(uy_array_wideextent)
            
            
          
      
        # Make plot
                        
        rms_total = np.asarray(rms_total)
        rms_total_extent = np.asarray(rms_total_extent)
        
        # Construct an image linearly increasing in y
        xv, yv = np.meshgrid(np.linspace(np.asarray(widths_new).min(),np.asarray(widths_new).max(),9), np.linspace(0,15,9))
        zv = yv
        
        z_new = abs(np.array([rms_total_extent-rms_total]))
        
        x_dist, y_dist = np.meshgrid(z_new, np.linspace(0,15,9))
        
        # Draw the image over the whole plot area
        plt.rcParams.update(params) 
        
        levels = np.linspace(x_dist.min(), x_dist.max(), 100)
        cs = ax1.contourf(xv,yv,x_dist,levels = levels,cmap='RdYlGn_r', extend = 'both', pad = 20, rotation = 180)
        
        # Erase above the data by filling with white
        ax1.fill_between(widths_new, rms_total, rms_total.min(), color='w')
        ax1.fill_between(widths_new, rms_total_extent, rms_total.max(), color='w')
        
        # Make the line plot over the top
        ax1.plot(widths_new, rms_total, 'k-', linewidth=2)
        ax1.plot(widths_new, rms_total_extent,'k--',linewidth=2)
        
        ax1.set_ylim(rms_total.min(), rms_total.max())
        
        cax = fig.add_axes([0.20,0.15,0.5,0.02])
        cbar = plt.colorbar(cs,cax=cax,orientation = 'horizontal')
        #cbar.set_label('distance') 
        cbar.set_ticks([50,100])
        cbar.set_label('Distance', labelpad = 4, fontsize = 10)
        cax.xaxis.set_label_position('top')
        
        ax1.set_xlabel('Channel widths' + ' [m]',labelpad=15)
        ax1.set_ylabel('RMS of hydrostatic deviation [m]',labelpad=15)
       
     
        ax1.spines['top'].set_visible(True)
        ax1.spines['bottom'].set_visible(True)
        ax1.spines['right'].set_visible(True)
        ax1.spines['left'].set_visible(True)
        ax1.tick_params(direction='in',length=6,width=2)


        # Make twinx velocity plot
        
           # Plt properties for ax2 (bridging)  
        ax2 = ax1.twinx()
        
        ax2.plot(widths_new, rms_uy_normal, 'r-')
        ax2.plot(widths_new, rms_uy_extent, 'r--')
        ax2.plot(widths_new, rms_uy_wideextent, 'r:')
        
        
        ax2.set_ylabel('Channel - velocity $u_{y}$ [m/a]',labelpad = 15)
        #ax2.set_ylim(0.3,0.32)  
        ax2.tick_params('y')
        
        
        legend_ax2 = ax2.legend(['Velocity(regular domain)','Velocity (extended domain)', 'Velocity (wide extended domain)'],loc=1)
        frame_ax2 = legend_ax2.get_frame()
        frame_ax2.set_facecolor('0.7')
        frame_ax2.set_edgecolor('0.7')
        

        legend = ax1.legend(['t = ' + str(t*5)+' a'+ ' regular','t = ' + str(t*5)+' a' + ' extent'],loc=2)
                
        frame = legend.get_frame()
        frame.set_facecolor('0.7')
        frame.set_edgecolor('0.7')
        ax2.tick_params(direction='in',length=6,width=2)
        

       
        path = str('plots/')
        fname= str('allwidths_corr_2d_alldomains'+ '_' + str(t*5) + '.eps')
        
        fig.savefig(path + fname, format = 'eps', dpi=1000)
            
            
            
          
    
    return rms_total, rms_total_extent, widths_new