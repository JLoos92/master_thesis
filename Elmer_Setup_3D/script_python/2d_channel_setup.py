#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:47:52 2019

@author: jloos
"""


from __main__ import ModelRun
from matplotlib.mlab import griddata

channel = ModelRun(900000,0,0,0,10,"vtu")

npoints = channel.npoints
points = channel.Points
dict_var_names = channel.dict_var_names
dict_var_values = channel.dict_var_values

x = points[:,0]
y = points[:,1]

sxy = channel.sxy
fs_lower = channel.fs_lower
fs_upper = channel.fs_upper 
#plt.plot(points[:,0],points[:,1],sxy, 'ro')
pointsize = 10

plt.scatter(x,y,pointsize,sxy)



fig = plt.figure(figsize = (20,15))
plt.tripcolor(x,y,sxy)



surface = channel.compute_hydrostatic_thickness()
new_calc = surface[1]

lower = surface [4]
upper = surface[3]

lower_2 = surface[5]

hh = surface[2]
calcccc = surface[1] 
caldddd = abs(surface[0])
calc = surface[1]
fig = plt.figure(figsize = (20,15))
plt.plot(calc)
plt.plot(calcccc)



p_w = 1000.0 # kg m−3 )
p_i = 910.0  

together = upper + abs(lower)

hallo = np.divide((p_w*upper),(p_w-p_i))

fig = plt.figure(figsize = (20,15))
plt.plot(hallo)
plt.plot(together)