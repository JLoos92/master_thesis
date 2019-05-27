#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:21:19 2019

@author: jloos
"""
from __main__ import ModelRun
from matplotlib.mlab import griddata
import numpy.ma as ma

fig = plt.figure(figsize = (20,15))

for t in range(49,50):
      
    channel = ModelRun(500000,0,0,0,t,"vtu")
    
    ht = channel.compute_hydrostatic_thickness()
    upper=ht[2]
    lower=ht[3]
    
    x = ht[0]
    
    calc = ht[1]
    calc = calc + upper
   
    plt.plot(x,lower,'r-')
    plt.plot(x,calc,'b-')
    
plt.show()

  


