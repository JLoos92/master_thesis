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
                 **kwargs): 
        
       # self.t = t
        
        list_widths = ModelRun(200000,0,0,0,20,"vtu").list_widths
        list_widths.sort(key=int)
        list_widths = list_widths[20:25]
        
        fig1 = plt.figure(figsize = (15,15))
        
        rms_total = []
        for i in range (ModelRun(20000,0,0,0,20,"vtu").num_timesteps):
            mr = ModelRun(50000,0,0,0,i,"vtu")
            
            ht = mr.compute_hydrostatic_thickness()
            
            
            # compute deviation
            self.lower = ht[3]
            self.calc_thickness_bs = ht[1]
            self.new_x = ht[0]
            self.hydrostatic_deviation = self.calc_thickness_bs - self.lower
            rms = np.sqrt(np.mean(self.hydrostatic_deviation**2))
                   
            rms_total.append(rms)
           
             
            plt.plot(rms_total)
            
            plt.xlabel('Years * 5')
            plt.ylabel('RMS of hydrostatic deviation [m]')
           
           
        fig1.show()