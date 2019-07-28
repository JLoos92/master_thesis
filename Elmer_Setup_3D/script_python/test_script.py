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













