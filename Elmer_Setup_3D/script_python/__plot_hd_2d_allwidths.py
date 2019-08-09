#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:33:32 2019

@author: jloos
"""

import pandas as pd
from main import ModelRun
from __plot_params import params_vertical

#
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter






def compute_allwidths_hd_2d(t1=None,
                            t2=None,
                            t3=None,
                            t4=None,
                        **kwargs):

    '''
    Computes rms of the hydrostatic deviation. The hydrostatic de-
    viation is saved in an 1d array for the 2d simulation case.       
    '''    
    
    
    # List of widths for extended and tight domain
    list_widths = list(np.arange(5000,10000,1000))
    list_widths_1 = list(np.arange(10000,100000,10000))
    list_widths_2 = list(np.arange(200000,900000,100000))
    
    list_widths_1.extend(list_widths_2)
    
    
    
    num_timesteps = ModelRun(150,20000,0,'extent',50,"2").num_timesteps
    num_timesteps = num_timesteps -2
   
    
    nums = ["(a)","(b)","(c)","(d)"]
    times = [t1, t2,t3,t4]
    # place text box in upper left in axes coords      
    props = dict(boxstyle='round', facecolor='wheat')
    
    # Define rc_params for figure
    nrow = 1; ncol = 4;
    fig, ax = plt.subplots(nrows=nrow, ncols=ncol,sharey = True, squeeze = False, constrained_layout=True)
    
    plt.subplots_adjust(bottom=0.1,wspace=0)
    
    
    for ax1,num,t in zip(ax.reshape(-1),nums,times): 
        
        
        # Hydrostatic deviation
        rms_total = []
        rms_total_extent = []
        rms_total_wideextent = []
        
        # Velocities
        rms_uy_normal = []
        rms_uy_extent = []
        rms_uy_wideextent = []
        
        peak = []
        peak_extent = []
        peak_wideextent = []
        
        time = []
        widths_new = []
    
    
        for width in list_widths_1:
            
            
            # Define width in m
            original_width = int(round(np.sqrt(width)*2))
            widths_new.append(original_width)
            
            # Set up class
            mr = ModelRun(150,width,0,0,t,"2")
            mr_extent = ModelRun(150,width,0,'extent',t,'2')
            mr_wideextent = ModelRun(150,width,0,'wideextent',t,'2')
            
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
            ht_wideextent = mr_wideextent.compute_hydrostatic_thickness()
            
            
            
            # compute deviation regular
            lower = ht[3]
            calc_thickness_bs = ht[1]
            
            #extent
            lower_extent = ht_extent[3]
            calc_thickness_bs_extent = ht_extent[1]
            
            #wideextent
            lower_wideextent = ht_wideextent[3]
            calc_thickness_bs_wideextent = ht_wideextent[1]

            
            hydrostatic_deviation = calc_thickness_bs - lower
            hydrostatic_deviation_extent = calc_thickness_bs_extent - lower_extent
            hydrostatic_deviation_wideextent = calc_thickness_bs_wideextent - lower_wideextent
        
            
            # Adapt ratio of wider domain
            x_extent = ht_extent[0]
            df_extent = pd.DataFrame({'hd_wide':hydrostatic_deviation_extent,'x_extent':x_extent})
            df_extent = df_extent[df_extent.x_extent>-10000]
            df_extent = df_extent[df_extent.x_extent<10000]
            hd_extent = df_extent.iloc[:,0].values
            
            # Adapt ratio of wider domain
            x_wideextent = ht_wideextent[0]
            df_wideextent = pd.DataFrame({'hd_wide':hydrostatic_deviation_wideextent,'x_wideextent':x_wideextent})
            df_wideextent = df_wideextent[df_wideextent.x_wideextent>-10000]
            df_wideextent = df_wideextent[df_wideextent.x_wideextent<10000]
            hd_wideextent = df_wideextent.iloc[:,0].values
            
            
            
            # Adapt ratio of extra-wide domain
        
            rms = np.sqrt(np.mean(hydrostatic_deviation**2))
            rms_extent = np.sqrt(np.mean(hd_extent**2))
            rms_wideextent = np.sqrt(np.mean(hd_wideextent**2))
            
           
            # Append rms value for timestep i
            rms_total.append(rms)
            rms_total_extent.append(rms_extent)
            rms_total_wideextent.append(rms_wideextent)
            
            # Append uy values
            rms_uy_normal.append(uy_array_normal)
            rms_uy_extent.append(uy_array_extent)
            rms_uy_wideextent.append(uy_array_wideextent)
            
            # Append peak deviation
            peak.append(ht[8])
            peak_extent.append(ht_extent[8])
            peak_wideextent.append(ht_wideextent[8])
            
            
          
      
        # Make plot
                        
        rms_total = np.asarray(rms_total)
        rms_total_extent = np.asarray(rms_total_extent)
        
        # Construct an image linearly increasing in y
        xv, yv = np.meshgrid(np.linspace(np.asarray(widths_new).min(),np.asarray(widths_new).max(),16), np.linspace(0,15,16))
        
        z_new = abs(np.array([rms_total_extent-rms_total]))
        
        x_dist, y_dist = np.meshgrid(z_new, np.linspace(0,15,16))
        
        # Draw the image over the whole plot area
        plt.rcParams.update(params_vertical) 
        
        levels = np.linspace(x_dist.min(), x_dist.max(), 100)
        #cs = ax1.contourf(xv,yv,x_dist,levels = levels,cmap='Blues_r', extend = 'both')
        
        # Erase above the data by filling with white
        ax1.fill_between(widths_new, rms_total, rms_total.min(), color='w')
        ax1.fill_between(widths_new, rms_total_extent, rms_total.max(), color='w')
        
        # Make the line plot over the top
        ax1.plot(widths_new, rms_total, 'k-', linewidth=1.5)
        ax1.plot(widths_new, rms_total_extent,'k--',linewidth=1.5)
        ax1.plot(widths_new, rms_total_wideextent,'k:',linewidth=1.5)
        
        
        
        ax1.text(0.39, 0.95, num + '  t = ' + str(t*5) + 'a', transform=ax1.transAxes, 
                verticalalignment='top', bbox=props, weight='bold',fontsize=8.5)  
        ax1.set_ylim(0, 8)
        ax1.set_yticks([1,3,5,7])
        ax1.set_xticks([500,1500])


        
        plt.setp(ax1.get_xticklabels(),fontweight = 'bold')
        plt.setp(ax1.get_yticklabels(),fontweight = 'bold')
     
        ax1.spines['top'].set_visible(True)
        ax1.spines['bottom'].set_visible(True)
        ax1.spines['right'].set_visible(True)
        ax1.spines['left'].set_visible(True)
        ax1.tick_params(direction='in',length=6,width=2)


        # Make twinx velocity plot
        
           # Plt properties for ax2 (bridging)  
        ax2 = ax1.twinx()
        
        ax2.plot(widths_new, peak, 'r-')
        ax2.plot(widths_new, peak_extent, 'r--')
        ax2.plot(widths_new, peak_wideextent, 'r:')
        
        # Set label and color for ax2 (second y-axis)
        
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
            
        ax2.tick_params('y')
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        
        
        
        ax2.tick_params(direction='in',length=6,width=2)
        
        ax2.set_yticks([1,10,20,30])
        ax2.set_ylim(0,35)
        ax2.set_yticklabels([])
       
        fig.add_subplot(ax1)
        
        
        
        
    ax[0][0].set_ylabel('RMS of hydrostatic deviation [m]') 
    legend = ax[0][3].legend(['regular domain','extended "" ', 'wide-extended "" '],loc="lower center",bbox_to_anchor=[0.52,0.65])                
    frame = legend.get_frame()
    frame.set_facecolor('0.7')
    frame.set_edgecolor('0.7')
    ax2.axes.get_yaxis().set_visible(True)
    ax2.set_yticks([0,10,20,30])
    ax2.set_yticklabels([0,10,20,30])
    ax2.set_xticks([500,1500])
    ax2.set_yticks([0,10,20,30])
    ax2.set_ylabel('Channel peak deviation [m]',visible = True, color= 'r')
    plt.setp(ax2.get_yticklabels(),fontweight = 'bold',color='r')
    #fig.suptitle('Widths vs. deviation @ multiple domains' + ' ', weight='bold', fontsize = 12,y=0.94)    
    fig.text(0.39,0,'Channel widths [m]',va = 'center',fontsize=11)
     
        
        
        
        
    path = str('plots/allwidths_2d_hd_velo/')
    fname_eps= str('TESTallwidths_corr_2d_alldomains'+ '_' + str(t*5) + '.pdf')
    fname_png= str('TESTallwidths_corr_2d_alldomains'+ '_' + str(t*5) + '.png')
        
    fig.savefig(path + fname_eps, format = 'pdf', dpi=1000,bbox_inches='tight')
    fig.savefig(path + fname_png, format = 'png', dpi=1000,bbox_inches='tight')
            
            
            
          
    
    return rms_total, rms_total_extent, widths_new
    