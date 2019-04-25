#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 14:43:45 2019

@author: jloos
"""


from pylab import rcParams
import numpy as np
import matplotlib.pyplot as plt
from __main__ import ModelRun
import pandas as pd
import math






def compute_correlation(t=None, 
                        **kwargs):


    #====================================================================== 
    # Mean of hydrostatic deviation ()
    #====================================================================== 
    
    
    # Defines region of interest
    x1 = 1060000
    x2 = 1080000
    y1 = -100 
    y2 = 100
    
    
    # Width of bedrock bumps: x2 (gauss-function)
    runs = [50,75,100,125,150,200,225,250,275,300,400,450,475,500]
    runs = np.asarray(runs)
    real_width = runs * 2
    
    hd_list_mean = []
    hd_list_rms = []
    hd_list_rms_noise = []
    hd_list_rms_noise_2 = []
    
    
    if t is None:
        t=200
    elif t is not None:
        t = t
    
    for i in runs:
    
        hydrostatic_thick = ModelRun(250,i,i,0,t).compute_hydrostatic_thickness()
        
        x = hydrostatic_thick[0]
        y = hydrostatic_thick[1]
        ht = hydrostatic_thick[2]
        
        df = pd.DataFrame({'X':x ,'Y':y, 'Hydrostatic deviation':ht})
        df_sorted = df.sort_values(by = ['Y'])
        
        df_sorted_new = df_sorted[(df_sorted['Y']>= y1)
                                & (df_sorted['Y']<= y2)]
        df_sorted_new = df_sorted_new[(df_sorted['X']>= x1)
                                & (df_sorted_new['X']<= x2)]
        
        hd_range = df_sorted_new.iloc[:,2].values
        
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
        
        
        # rms method for calculating noise
        def rmsvalue(arr,n):
            square = 0
            mean = 0.0
            root= 0.0
            
            for i in range(0,n):
                square += (arr[i]**2)
                
                mean = (square / (float)(n))
                
                root = math.sqrt(mean)
                
                return root
        
        n = len(hd_range)
        hd_rms = rmsvalue(hd_range,n)
        hd_list_rms.append(hd_rms)
        
        n = len(hd_range_noise)
        hd_rms_noise = rmsvalue(hd_range_noise,n)
        hd_list_rms_noise.append(hd_rms_noise)
        
        n = len(hd_range_noise_2)
        hd_rms_noise_2 = rmsvalue(hd_range_noise_2,n)
        hd_list_rms_noise_2.append(hd_rms_noise_2)
    
    
        
        
    hd_array_mean = np.asarray(hd_list_mean)    
    hd_array_rms = np.asarray(hd_list_rms) 
    hd_array_rms_noise = np.asarray(hd_list_rms_noise) 
    hd_array_rms_noise_2 = np.asarray(hd_list_rms_noise_2) 
        
    hd_new = hd_array_rms - (hd_array_rms_noise+hd_array_rms_noise_2)
    
    
    plt.plot(real_width,hd_array_mean)
    plt.plot(real_width,hd_array_mean,'bo')
    plt.xlabel(' Width of bedrock bumps [m] ')
    plt.ylabel('Mean hydrostatic deviation in boundaries')
    plt.title('Hydrostatic deviation as a function of bedrock-bump width')
    plt.show()
    
    plt.plot(real_width,hd_array_rms,'bo')
    plt.xlabel('Width of bedrock bumps [m] ')
    plt.ylabel('RMS hydrostatic deviation in boundaries')
    plt.title('Width')
    plt.show()
    
    
    plt.plot(real_width,hd_new,'bo')
    plt.xlabel(' Width of bedrock bumps [m] ')
    plt.ylabel('RMS hydrostatic deviation in boundaries')
    plt.title('Without noise')
    plt.show()
    
   





        