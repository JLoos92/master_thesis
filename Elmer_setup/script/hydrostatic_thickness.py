#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 16:21:44 2019

@author: Schmulius
"""
import numpy as np

from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
import plotly as py
import turbulucid


p_w = 1027 #kg m−3 ), ice (ρi =918kgm−3), and air (ρa =2kgm−3):
p_i = 900
p_a = 2
H_a = 0


Points_0 =  vtk_to_numpy(xmlReader2.GetOutput().GetPoints().GetData())
Points_60 =  vtk_to_numpy(xmlReader60.GetOutput().GetPoints().GetData())
velo =  vtk_to_numpy(xmlReader2.GetOutput().GetPointData().GetArray(17))

H2 = ((p_w*Points_0[:,2])/(p_w-p_i))-(((p_a-p_i)/(p_i-p_w))* H_a)
H60 = ((p_w*Points_60[:,2])/(p_w-p_i))-(((p_a-p_i)/(p_i-p_w))*H_a)

H_thickness = H2 - H60
H_thickness_new = H2 - Points_60[:,2]

#fig = plt.imshow(H_thickness, extent=[min(Points_0[:,0]),max(Points_0[:,0]),min(Points_0[:,1]),max(Points_0[:,1])],origin="lower")

#imagesc(Points_0[:,0], Points_0[:,1], H_thickness)
x = Points_0[:,0]
y = Points_0[:,1]
z = Points_0[:,2]
velo = velo[:,0]
x = np.expand_dims(x, axis=0)
y = np.expand_dims(y, axis=0)
z = np.expand_dims(z, axis=0)
H_thickness = np.expand_dims(H_thickness, axis=0)
velo = np.expand_dims(velo, axis=0)


#fig = plt.imshow(velo, extent=[min(Points_0[:,0]),max(Points_0[:,0]),min(Points_0[:,1]),max(Points_0[:,1])],origin="lower",interpolation='bessel')
#plt.colorbar()
