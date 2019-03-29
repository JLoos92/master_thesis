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
#import plotly as py
import turbulucid
from modelrun import *
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
import numpy as np
from tempfile import TemporaryFile
import numpy as np
import matplotlib
from pylab import *
import numpy as np
from pandas import DataFrame, Series
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
from scipy.interpolate import griddata,interp2d
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay,ConvexHull
from matplotlib.tri import Triangulation,TriInterpolator, LinearTriInterpolator,CubicTriInterpolator, TriFinder
import pylab
from numpy.linalg import inv
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import LinearNDInterpolator,Rbf, UnivariateSpline,CloughTocher2DInterpolator
import concave_hull


# Setup for figure size
matplotlib.get_backend()
matplotlib.rcParams['figure.figsize'] = (25,15)




# Use class and method to compute hydrostatic thickness at timestep

h_thickness0_250 = ModelRun(250,250,250,0,0).compute_hydrostatic_thickness()
h_thickness25 = ModelRun(250,250,250,0,25).compute_hydrostatic_thickness()
h_thickness50_250 = ModelRun(250,250,250,0,50).compute_hydrostatic_thickness()
h_thickness75 = ModelRun(250,250,250,0,75).compute_hydrostatic_thickness()
h_thickness100 = ModelRun(250,250,250,0,100).compute_hydrostatic_thickness()

#h_thickness0_500 = ModelRun(250,500,500,0,0).compute_hydrostatic_thickness()
h_thickness50_500 = ModelRun(250,500,500,0,50).compute_hydrostatic_thickness()

#ht_0_250 = h_thickness0_250[2] 
ht_50_250= h_thickness50_250[2]

#ht_0_500 = h_thickness0_500[2] 
ht_50_500= h_thickness50_500[2]


x = h_thickness50_250[0]
y = h_thickness50_250[1]

points = [x,y]
points = np.asarray(points)
points = points.transpose()
#cutter = ModelRun(250,250,250,0,50).cutter()
#cutter100 = ModelRun(250,250,250,0,100).cutter()
#
#
#points = vtk_to_numpy(cutter.GetOutput().GetPoints().GetData())
##_cutter = vtk_to_numpy(cutter.GetPoints().GetData())
#zs = vtk_to_numpy(cutter.GetOutput().GetPointData().GetArray('zs'))
#zb = vtk_to_numpy(cutter.GetOutput().GetPointData().GetArray('zb'))
#
#zs100 = vtk_to_numpy(cutter100.GetOutput().GetPointData().GetArray('zs'))
#zb100 = vtk_to_numpy(cutter100.GetOutput().GetPointData().GetArray('zb'))
#zeros_zs = np.nonzero(zs100)
#zeros_zb = np.nonzero(zb100)



# Triangulation or delaunay
tri = Triangulation(x,y)
delaunay = Delaunay(points)

# Plot mesh

plt.triplot(x,y,delaunay.simplices,alpha=0.2)
plt.tripcolor(tri,ht_50_500,shading='gouraud')
plt.title('Delaunay triangulation of refined mesh')
plt.clim(-5,10)
plt.colorbar()
#plt.gca().set_aspect('equal')
plt.show()

plt.triplot(x,y,delaunay.simplices,alpha=0.2)
plt.tripcolor(tri,ht_50_250,shading='gouraud')
plt.title('Delaunay triangulation of refined mesh')
plt.clim(-5,10)
plt.colorbar()
#plt.gca().set_aspect('equal')
plt.show()