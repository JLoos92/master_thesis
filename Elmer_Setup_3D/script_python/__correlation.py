#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 14:43:45 2019

@author: jloos
"""



import numpy as np
import matplotlib.pyplot as plt
from __main__ import ModelRun
import pandas as pd
import math
from __plot_params import params_3d,params_horizontal

import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'





def compute_correlation_3d(t=None, 
                        **kwargs):

    ''''
    
    
    '''

    #====================================================================== 
    # Mean of hydrostatic deviation ()
    #====================================================================== 
    
    
    # Defines region of interest
    x1 = 1060000
    x2 = 1080000
    y1 = -5000 
    y2 = 5000
    
    
    # Width of bedrock bumps: x2 (gauss-function)
    runs = [100,125,150,200,225,250,275,300,400,450,475,500,600]
    
    runs = np.asarray(runs)
    real_width = runs * 2
    
    

    if t is None:
        t=200
    elif t is not None:
        t = t
    # Time + numbers for ax
    nums = ["(a)","(b)","(c)","(d)","(e)","(f)"]
    
    orange = '#D55E00'
    fig, axs = plt.subplots(1,4, sharex=True, sharey = True) 
    times = [30,70,100,200]
    
     # Custom model load from __plot_params
    plt.rcParams.update(params_horizontal)
    plt.subplots_adjust(hspace=0,wspace=0)
    
    for ax1,num,t in zip(axs.flat,nums,times):
        hd_list_mean = []
        hd_list_rms = []
        hd_list_rms_noise = []
        hd_list_rms_noise_2 = []
        for i in runs:
        
            hydrostatic_thick = ModelRun(250,i,i,0,t).compute_hydrostatic_thickness()
            
            x = hydrostatic_thick[0]
            y = hydrostatic_thick[1]
            ht = hydrostatic_thick[2]
            
            df = pd.DataFrame({'X':x ,'Y':y, 'Hydrostatic deviation':ht})
            df_sorted = df.sort_values(by = ['Y'])
            
            df_sorted_new = df_sorted[(df_sorted['Y']>=  y1)
                                    & (df_sorted['Y']<= y2)]
            df_sorted_new = df_sorted_new[(df_sorted['X']>= x1)
                                    & (df_sorted_new['X']<= x2)]
            
            hd_range = df_sorted_new.iloc[:,2].values
            #print('hd_range',hd_range)
            hd_mean = np.mean(hd_range)
            
            hd_list_mean.append(hd_mean)
        
        
        #====================================================================== 
        # Noise of area
        #====================================================================== 
        
            df_sorted_noise = df_sorted[(df_sorted['Y']>= y2)]
            df_sorted_noise = df_sorted_noise[(df_sorted_noise['X']>= x1)
                                            & (df_sorted_noise['X']<= x2)]
            hd_range_noise = df_sorted_noise.iloc[:,2].values
            
        
            df_sorted_noise_2 = df_sorted[(df_sorted['Y']<= y1)]
            df_sorted_noise_2 = df_sorted_noise_2[(df_sorted_noise_2['X']>= x1)
                                                & (df_sorted_noise_2['X']<= x2)]
            
            hd_range_noise_2 = df_sorted_noise_2.iloc[:,2].values
            
            
            # rms 
            def rmsvalue(arr,n):
                square = 0
                mean = 0.0
                root= 0.0
                
                for i in range(0,n):
                    square += (arr[i]**2)
                    #print(square)
                    mean = (square / (float)(n))
                    #print('mean = ', mean)
                    root = math.sqrt(mean)
                    
                    return root
                
            rms = np.sqrt(np.mean(hd_range**2))
            n = len(hd_range)
            hd_rms = rmsvalue(hd_range,n)
            hd_list_rms.append(rms)
            print('rms = ', rms)
            
            rms_noise = np.sqrt(np.mean(hd_range_noise**2))
            n = len(hd_range_noise)
            hd_rms_noise = rmsvalue(hd_range_noise,n)
            hd_list_rms_noise.append(rms_noise)
           
            
            rms_noise_2 = np.sqrt(np.mean(hd_range_noise_2**2))
            n = len(hd_range_noise_2)
            hd_rms_noise_2 = rmsvalue(hd_range_noise_2,n)
            hd_list_rms_noise_2.append(rms_noise_2)
            
        
            
            
        hd_array_mean = np.asarray(hd_list_mean)    
        hd_array_rms = np.asarray(hd_list_rms)
        
        hd_array_rms_noise = np.asarray(hd_list_rms_noise) 
        hd_array_rms_noise_2 = np.asarray(hd_list_rms_noise_2) 
            
        hd_new = hd_array_rms - (hd_array_rms_noise+hd_array_rms_noise_2)
        hd_noise_ensemble = (hd_array_rms_noise+hd_array_rms_noise_2)
        
         # place text box in upper left in axes coords      
        props = dict(boxstyle='round', facecolor='wheat')
        

        ax1.text(0.2, 0.9,str(num) + 't = ' + str(t) + 'a', transform=ax1.transAxes,
                 verticalalignment='top', bbox=props,weight='bold')
        
        ax1.plot(real_width,hd_array_rms,'k-')
        ax1.grid(True)
        
        fig.add_subplot(ax1)
    
    fig.text(0.4,0.05,'Channel widths [m]',va = 'center',fontsize=6.5)
    fig.text(0.03,0.5,'RMS of hydrostatic deviation [m]',va = 'center',rotation = 'vertical',fontsize=6.5)
     # Save figures    
    path = str('plots/Final_plots/')
       
    fname_png = str('correlation_3d_' + str(t) + '.png')
    fname_pdf = str('correlation_3d_' + str(t) + '.pdf')

    
    plt.savefig(path + fname_png, format = 'png',dpi=1000,bbox_inches='tight')
    plt.savefig(path + fname_pdf, format = 'pdf',dpi=1000,bbox_inches='tight')  
    
    
    
    plt.show()
    
    return hd_array_rms



    

    
    