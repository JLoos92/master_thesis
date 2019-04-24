#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 15:31:28 2019

@author: Schmulius
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#matplotlib.use('nbagg')
from pylab import *
import numpy as np
from pandas import DataFrame, Series
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
#import plotly as py
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

p_w = 1027 #kg m−3 ), ice (ρi =918kgm−3), and air (ρa =2kgm−3):
p_i = 900
p_a = 2
H_a = 0


Points_0 =  vtk_to_numpy(case_0.points)
#Points_60 =  vtk_to_numpy(xmlReader60.GetOutput().GetPoints().GetData())


#H2 = ((p_w*Points_0[:,2])/(p_w-p_i))-(((p_a-p_i)/(p_i-p_w))* H_a)
#H60 = ((p_w*Points_60[:,2])/(p_w-p_i))-(((p_a-p_i)/(p_i-p_w))*H_a)

#H_thickness = H2 - H60
#H_thickness_new = H2 - Points_60[:,2]


ly_max = 5000
ly_min = -5000
lx_min= 1010000
lx_max= 1080000

xv = np.linspace(lx_min, lx_max, num=1001)
yv = np.linspace(ly_min, ly_max, num=201)

velo =  vtk_to_numpy(xmlReader2.GetOutput().GetPointData().GetArray(17))
dsdt60 =  vtk_to_numpy(xmlReader60.GetOutput().GetPointData().GetArray(2))
dsdt2 =  vtk_to_numpy(xmlReader2.GetOutput().GetPointData().GetArray(2))
stress60 =  vtk_to_numpy(xmlReader60.GetOutput().GetPointData().GetArray(18))

yy,xx = np.meshgrid(yv,xv)

velo = velocity[:,1]
stress_zz = stress60[:,2]
stress_xx = stress60[:,0]
stress_yy = stress60[:,1]

llx = 100000
lly = 20000

ZS = genfromtxt('ZS.xyz')
Points2 = np.load('Points2.npy')

x = Points2[:,0]
y = Points2[:,1]
z = Points2[:,2]
#x = arange(0,10)
#y = exp(-x/3.0)
points = np.delete(Points2, 2, 1) 


grid_z0 = griddata(points, stress_xx, (xx, yy), method='linear')
grid_z1 = griddata(points, stress_zz, (xx, yy), method='linear')
grid_z2 = griddata(points, H_thickness, (xx, yy), method='linear')
grid_z3 = griddata(points, z, (xx, yy), method='linear')




plt.subplot(221)
plt.imshow(grid_z0, cmap='jet')
#pcm = ax[0].pcolor(X, Y, Z1, norm=colors.LogNorm(vmin=Z1.min(), vmax=Z1.max()),cmap='jet')
#fig.colorbar(pcm, ax=ax[0], extend='max')

plt.title('stress_xx tensor')
#plt.colorbar()


plt.subplot(222)
plt.imshow(grid_z1, cmap='jet')
#pcm = ax[0].pcolor(X, Y, Z1, norm=colors.LogNorm(vmin=Z1.min(), vmax=Z1.max()),cmap='jet')
#fig.colorbar(pcm, ax=ax[0], extend='max')

plt.title('stress_zz tensor')
plt.colorbar()
plt.subplot(223)
plt.imshow(grid_z2, cmap='jet')
#pcm = ax[0].pcolor(X, Y, Z1, norm=colors.LogNorm(vmin=Z1.min(), vmax=Z1.max()),cmap='jet')
#fig.colorbar(pcm, ax=ax[0], extend='max')
plt.title('Hydrost. thickness')
plt.colorbar()
plt.subplot(224)
plt.imshow(grid_z3, cmap='jet')
#pcm = ax[0].pcolor(X, Y, Z1, norm=colors.LogNorm(vmin=Z1.min(), vmax=Z1.max()),cmap='jet')
#fig.colorbar(pcm, ax=ax[0], extend='max')

plt.title('stress_yy tensor')

plt.show()