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
    y1 = -100 
    y2 = 100
    
    
    # Width of bedrock bumps: x2 (gauss-function)
    runs = [50,75,100,125,150,200,225,250,275,300,400,450,475,500,600]
    
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
        print(rms_noise)
        
        rms_noise_2 = np.sqrt(np.mean(hd_range_noise_2**2))
        n = len(hd_range_noise_2)
        hd_rms_noise_2 = rmsvalue(hd_range_noise_2,n)
        hd_list_rms_noise_2.append(rms_noise_2)
        print(rms_noise_2)
    
        
        
    hd_array_mean = np.asarray(hd_list_mean)    
    hd_array_rms = np.asarray(hd_list_rms)
    print(hd_array_rms)
    hd_array_rms_noise = np.asarray(hd_list_rms_noise) 
    hd_array_rms_noise_2 = np.asarray(hd_list_rms_noise_2) 
        
    hd_new = hd_array_rms - (hd_array_rms_noise+hd_array_rms_noise_2)
    hd_noise_ensemble = (hd_array_rms_noise+hd_array_rms_noise_2)
    #print(hd_array_rms)
    print('noise_ensemble = ', hd_noise_ensemble)
    print('real area = ', hd_new)
    
#    plt.plot(real_width,hd_array_mean)
#    plt.plot(real_width,hd_array_mean,'bo')
#    plt.xlabel(' Width of bedrock bumps [m] ')
#    plt.ylabel('Mean hydrostatic deviation in boundaries')
#    plt.title('Hydrostatic deviation as a function of bedrock-bump width')
#    plt.show()
#    
#    plt.plot(real_width,hd_array_rms,'bo')
#    plt.xlabel('Width of bedrock bumps [m] ')
#    plt.ylabel('RMS hydrostatic deviation in boundaries')
#    plt.title('Width')
#    plt.show()
    
    
    
    
    
    fig = plt.figure(figsize = (15,10)) 
    plt.plot(real_width,hd_array_rms,'bo')
    plt.xlabel('Width of bedrock bumps [m] ',fontdict=font_label)
    plt.ylabel('RMS: Hydrostatic deviation [m]',fontdict=font_label)
    plt.title('',fontdict=font_title)
    plt.grid(True)
    plt.show()
    
   


def compute_correlation_2d(t=None, 
                        **kwargs):

    '''
    Computes root mean square of the hydrostatic deviation. The hydrostatic de-
    viation is saved in an 1d array for the 2d simulation case.       
    '''    
    
 
    
    list_widths = ModelRun(200000,0,0,0,50,"vtu").list_widths
    list_widths.sort(key=int)
    list_widths = list_widths[0:10]
    num_timesteps = ModelRun(20000,0,0,0,50,"vtu").num_timesteps
    
    
    fig1 = plt.figure(figsize = (15,15))
    
    for widths in list_widths:
        rms_total = []
        for i in range(1,100):
            
            mr = ModelRun(widths,0,0,0,i,"vtu")
            
            ht = mr.compute_hydrostatic_thickness()
            
            
            # compute deviation
            lower = ht[3]
            calc_thickness_bs = ht[1]
            
            hydrostatic_deviation = calc_thickness_bs - lower        
            rms = np.sqrt(np.mean(hydrostatic_deviation**2))
                   
            rms_total.append(rms)
            
            plt.plot(rms_total,label = 'Channel width =' + widths)
   
    
    
    plt.xlabel('Timesteps')
    plt.ylabel('RMS of hydrostatic deviation [m]')
    
    
    return rms_total





    
def plot_correlation_widths_2d(amp,
                               t,                               
                        **kwargs):

    '''
    Computes root mean square of the hydrostatic deviation. The hydrostatic de-
    viation is saved in an 1d array for the 2d simulation case.       
    '''    
    
    
    #==========================================================================
    # Setup fonts
    #==========================================================================
        

    font_annotation = {'color':'black',
                   'size': '17'
                   }    
   
    font_title = {'color':'black',
           'size':'20',
           'weight':'bold'
           }

    font_axes = {'color':'black',
           'size':'14',
           'weight':'normal'
           }
    
    font_label = {'color':'black',
           'size':'15',
           'weight':'normal'
           }
 
    
    list_widths = ModelRun(amp,20000,0,0,t,"2").df_amps_widths
    list_widths= list_widths[list_widths['Amplitudes'].isin([str(amp)])]
    list_widths = list_widths['Widths'].tolist()
    
    list_widths.sort(key=int)
    list_widths = list_widths
   
    num_timesteps = ModelRun(amp,20000,0,0,t,"2").num_timesteps
    original_widths = []
    for i in range(len(list_widths)):
        original_widths_new = int(list_widths[i])
        original_widths.append(original_widths_new)
        
    original_widths = np.sqrt(original_widths)*2*2
    fig1 = plt.figure(figsize = (15,15))
    rms_total = []
    for width in list_widths:
        
            

            mr = ModelRun(amp,width,0,0,t,"2")
            
            ht = mr.compute_hydrostatic_thickness()
            
            
            # compute deviation
            lower = ht[3]
            calc_thickness_bs = ht[1]
            
            hydrostatic_deviation = calc_thickness_bs - lower        
            rms = np.sqrt(np.mean(hydrostatic_deviation**2))
                   
            rms_total.append(rms)
            
    plt.plot(original_widths,rms_total,'r*')
   
    plt.xlabel('Channel width [m]',fontdict = font_label)
    plt.ylabel('RMS of hydrostatic deviation [m]',fontdict = font_label)
    plt.legend(loc = 'upper left')
    
    
    