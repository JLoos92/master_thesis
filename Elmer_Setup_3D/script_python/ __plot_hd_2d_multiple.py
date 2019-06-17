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





class Plot_hydrostatic_deviation_2d_multiple():
    
    """
    Plots hydrostatic deviation as an 1d-array. Multiple = all runs of dictionary
    """
    

    
    
    
    def __init__(self,
                 width): 
        
       # self.t = t
        
        list_widths = ModelRun(150,20000,0,0,20,"2").list_widths
       
        
        fig1 = plt.figure(figsize = (15,15))
        
        rms_total = []
        rms_total_ex = []
        
        
        
        for i in range (ModelRun(150,20000,0,0,20,"2").num_timesteps):
            mr = ModelRun(150,width,0,0,i,"2")
            mr_extent = ModelRun(150,90000,0,'extent',i,"2")
            
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
           
             
            plt.plot(rms_total,'b-')
            plt.plot(rms_total_ex,'r-')
            
            plt.xlabel('Years * 5')
            plt.ylabel('RMS of hydrostatic deviation [m]')
           
           
        fig1.show()