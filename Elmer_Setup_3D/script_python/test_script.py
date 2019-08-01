# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 15:55:13 2018
"""



import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from scipy.interpolate import griddata
import pandas as pd



points_extent =  ModelRun(150,20000,0,'extent',10,"2").Points

x_extent = points_extent[:,0]
uy_extent = ModelRun(150,20000,0,'extent',10,'2').get_scalar('velocity')  
uy_extent = uy_extent[:,1] 

# Extended
df_extent = pd.DataFrame({'x_extent':x_extent,'uy_extent':uy_extent})
df_extent = df_extent[df_extent.x_extent>0]
df_extent = df_extent[df_extent.x_extent<10000]
sxy_array_extent = df_extent.iloc[:,1].values
rms_uy_extent = np.sqrt(np.mean(sxy_array_extent**2))





hd_extended = ModelRun(150,90000,0,0,200,"2").compute_hydrostatic_thickness()
ht_calc = hd_extended[1]
ht_model = hd_extended[3]

ht_calc_ind = int(abs(ht_calc.size/2)) 
ht_model_ind = int(abs(ht_model.size/2))

peak_calc = abs(ht_calc[ht_calc_ind])
peak_model = abs(ht_model[ht_model_ind])

# absolute difference (peak difference in m)
abs_diff = peak_calc-peak_model



hd_extended = ModelRun(150,90000,0,'extent',200,"2")




hd_extended = ModelRun(150,90000,0,0,200,"2").compute_hydrostatic_thickness()

peak = hd_extended[7]
