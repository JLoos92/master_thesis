#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:29:50 2019

@author: jloos
"""
from scipy import interpolate
from scipy.interpolate import griddata

mr = ModelRun(200000,0,0,0,50,"vtu")
            
ht = mr.compute_hydrostatic_thickness()

points = ht[4]
x = points[:,0]
y = points[:,1]
lower = ht[3]
new_points = x,y
xline = ht[0]
xx,yy = np.meshgrid(x,y)

sxy = mr.sxy
grid = griddata(new_points,sxy,(xx,yy),method='nearest')

plt.pcolormesh(xx,yy,grid,cmap = 'RdBu')
plt.plot(xline,lower,'b-')
