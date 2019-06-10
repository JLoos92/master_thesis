#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:29:50 2019

@author: jloos
"""
from scipy import interpolate
from scipy.interpolate import griddata
import numpy.ma as ma

mr = ModelRun(200000,0,0,0,50,"2")
            
ht = mr.compute_hydrostatic_thickness()

points = ht[4]
x = points[:,0]
y = points[:,1]
new_points = x,y

lower = ht[3]

xline = ht[0]


full_lower = ht[5]
full_upper = ht[6]


masked_fulllower =  ma.masked_where(full_lower>=0,full_lower)
masked_total = ma.masked_where(masked_fulllower<y,y)


sxy = mr.sxy

xx,yy = np.meshgrid(x,y)
grid = griddata(new_points,sxy,(xx,yy),method='nearest')


plt.pcolormesh(xx,yy,grid,cmap = 'RdBu')
plt.plot(xline,lower,'b*')
