#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:28:51 2019

@author: jloos
"""

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
from pylab import rcParams
import numpy as np
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
import turbulucid
from __main__ import ModelRun
from scipy.spatial import Delaunay,ConvexHull
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
from math import trunc
import numpy.ma as ma
from __plot_params import params




class Plot_hydrostatic_deviation_2d_multiple():
    
    """
    Plots hydrostatic deviation as an 1d-array. Multiple = all runs of dictionary
    """
    

    
    
    
    def __init__(self,
                 width_1,
                 width_2): 
        
       # self.t = t
        
        list_widths = ModelRun(150,20000,0,0,20,"2").list_widths
       
        
       
        fig1, ax1 = plt.subplots()  
        time = []
        original_width_1 = int(round(np.sqrt(width_1)*2)*2)
        original_width_2 = int(round(np.sqrt(width_2)*2)*2)
        rms_total = []
        rms_total_ex = []
        
        
        
        for i in range (3,100):
            mr = ModelRun(150,width_1,0,0,i,"2")
            mr_extent = ModelRun(150,width_2,0,0,i,"2")
            
            ht = mr.compute_hydrostatic_thickness()
            ht_extent = mr_extent.compute_hydrostatic_thickness()
            
            
            # compute deviation extent
            self.lower_ex = ht_extent[3]
            self.calc_thickness_bs_ex = ht_extent[1]
            self.new_x_ex = ht_extent[0]
            self.hydrostatic_deviation_extent = self.calc_thickness_bs_ex - self.lower_ex
            rms_ex = np.sqrt(np.mean(self.hydrostatic_deviation_extent**2))
                   
            rms_total_ex.append(rms_ex)
            
          # compute deviation regular
            self.lower = ht[3]
            self.calc_thickness_bs = ht[1]
            self.new_x = ht[0]
            self.hydrostatic_deviation = self.calc_thickness_bs - self.lower
            rms = np.sqrt(np.mean(self.hydrostatic_deviation**2))
                   
            rms_total.append(rms)
            time.append(i*5)
           
            # Make plot
                           
            
            ax1.plot(time,rms_total,'g-',linewidth=2)
            ax1.plot(time,rms_total_ex,'b-',linewidth=2)
            legend = ax1.legend(['channel width = ' + str(original_width_1)+' m', 'channel width = ' + str(original_width_2)+' m'],loc=1)
                    
            
                    
            frame = legend.get_frame()
            frame.set_facecolor('0.7')
            frame.set_edgecolor('0.7')
            
            
            # Plt properties
            plt.rcParams.update(params)
            ax1.spines['top'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['left'].set_visible(False)
            ax1.tick_params(direction='in',length=6,width=2)
            
            #plt.title('RMS of hd with multiple channel widths', y = 1.05)
            ax1.set_xlim(0,500)
            ax1.set_ylim(0,10)       
            ax1.set_xlabel('Time [a]',labelpad=20)
            ax1.set_ylabel('RMS of hydrostatic deviation [m]',labelpad=20)
            ax1.grid(linestyle = '--')
            
            
            path = str('plots/')
            fname= str('corr_2d_tightvswide' + str(original_width_1)+ str(original_width_2) + '.eps')
            
            fig1.savefig(path + fname, format = 'eps', dpi=1000,bbox_inches='tight')
            
           
        fig1.show()