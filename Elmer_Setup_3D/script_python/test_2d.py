#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:21:19 2019

@author: jloos
"""
from __main__ import ModelRun
from matplotlib.mlab import griddata
import numpy.ma as ma


fig = plt.figure(figsize = (20,10))  
channel = ModelRun(9000000,0,0,0,1,"vtu")

ht = channel.compute_hydrostatic_thickness()
upper=ht[6]
lower=ht[5]
points=ht[4]
x = ht[3]
calc = ht[1]

plt.plot(x,upper)
plt.plot(x,lower)
plt.plot(x,calc)


p_w = 1000.0
p_i = 910.0
thick_calc = -1 *( np.divide((p_w*upper),(p_w-p_i)))  
plt.plot(x,thick_calc) 
  


