#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 14:33:54 2019

@author: jloos
"""
from main import ModelRun
from __plot_params import params_bridging_2d
from scipy.signal import find_peaks_cwt


from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
import matplotlib.font_manager
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy.interpolate import griddata,interpolate


def local_minimum_1d(width,time):
    
    mr = ModelRun(150,width,0,0,time,"2")
                       
    ht = mr.compute_hydrostatic_thickness()

    lower = ht[3]
    calc_thickness_bs = ht[1]


    
    hydrostatic_deviation = calc_thickness_bs - lower

    peak = np.where(calc_thickness_bs==calc_thickness_bs.min())

    hd_peak = hydrostatic_deviation[peak]
    
    
    return hd_peak
 







   
