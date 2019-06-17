#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:33:32 2019

@author: jloos
"""

from __main__ import ModelRun
import pandas as pd
from scipy.interpolate import griddata

points = ModelRun(150,20000,0,0,50,"2").Points
sxy =  ModelRun(150,20000,0,0,50,"2").sxy
x = points[:,0]
y = points[:,1]



df = pd.DataFrame({'sxy':sxy,'x':x})
df = df[df.x>0]
df = df[df.x>5000]
sxy_array = df.iloc[:,0].values
rms = np.sqrt(np.mean(sxy_array**2))



#
#def perc(data):
#   median = np.zeros(data.shape[0])
#   perc_25 = np.zeros(data.shape[0])
#   perc_75 = np.zeros(data.shape[0])
#   for i in range(0, len(median)):
#       median[i] = np.median(data[:, i])
#       perc_25[i] = np.percentile(data[:, i], 25)
#       perc_75[i] = np.percentile(data[:, i], 75)
#   return median, perc_25, perc_75