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
channel = ModelRun(5000000,0,0,0,10,"vtu")

ht = channel.compute_hydrostatic_thickness()
upper=ht[6]
lower=ht[5]
points=ht[4]
x = ht[3]
plt.plot(x,upper)
plt.plot(x,lower)

sxy = channel.sxy
 
points = ht[4]
x = points[:,0]
y = points[:,1]
plt.tripcolor(x,y,sxy)
plt.show()

fig = plt.figure(figsize = (20,10)) 
indlow = ht[7]
indup = ht[8]

new = ma.masked_invalid(sxy)

plt.tripcolor(x,y,sxy)
plt.minorticks_on()
plt.show()


plt.plot(ht[5])
plt.plot(ht[1])
plt.plot(points[:,1])

