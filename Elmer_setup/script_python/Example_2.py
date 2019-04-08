#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 12:59:31 2018

@author: Schmulius
"""

import numpy as np
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
import turbulucid
from modelrun import ModelRun
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
from scipy.interpolate import griddata,interp2d,interp1d
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay,ConvexHull
from matplotlib.tri import Triangulation,TriInterpolator, LinearTriInterpolator,CubicTriInterpolator, TriFinder
import pylab
from numpy.linalg import inv
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import LinearNDInterpolator,Rbf, UnivariateSpline,CloughTocher2DInterpolator, griddata, CubicSpline
import concave_hull
from collections import OrderedDict
from matplotlib.colors import LogNorm
from pylab import rcParams


#hi = ModelRun(250,150,150,99,1).compute_concavehull(1069000,1068000,1070000,3)


#
#f_name = "/Volumes/esd01/docs/jloos/data_small/runs_elmerice_fixed/Mesh250_150150_0/Mesh/ForwardGLFixed0010.pvtu"
#
#xmlReader = vtk.vtkXMLPUnstructuredGridReader()
#xmlReader.SetFileName(f_name)
#xmlReader.Update()
#
#
#xcoord_1 = 1070000
#xcoord_2 = 1071000
#xcoord_3 = 1072000
#
#
#coordsx = [xcoord_1,xcoord_2,xcoord_3]
#
#for x in coordsx:
#    
#    line=vtk.vtkPlane()
#    lines = []
#    lines = line.SetOrigin(x,0,0)
#    line.SetNormal(1,0,0)
#                
#    
#    cutter = vtk.vtkCutter()   
#    cutter.SetCutFunction(line)
#    cutter.SetInputConnection(xmlReader.GetOutputPort())
#    cutter.Update()
#    line_points = vtk_to_numpy(xmlReader.GetOutput().GetPoints().GetData())
#                
#                # Rearrange matrix to fit input shape ()
#    hull_points = np.delete(line_points, 0, 1)
#                
#                # Compute concave hull at crosssection of line points
#    hull = concave_hull.compute(hull_points,3)
#    
#    
    
hydrostatic_thicknesses = {
    '1_ht_150_10' : ModelRun(250,150,150,'double',10).compute_hydrostatic_thickness(),
    '2_ht_150_25' : ModelRun(250,150,150,'double',25).compute_hydrostatic_thickness(),
    '3_ht_150_50' : ModelRun(250,150,150,'double',30).compute_hydrostatic_thickness(),
    '4_ht_150_100' : ModelRun(250,150,150,'double',40).compute_hydrostatic_thickness()
                            }
    
hydrostatic_thicknesses = OrderedDict(hydrostatic_thicknesses)    



# Define properties for subplots
    
font_title = {'color':'black',
       'size':'40',
       'weight':'bold'
            }

font_axes = {'color':'black',
       'size':'14',
       'weight':'italic'
           }

font_annotation = {'color':'black',
                   'size': '10'
                   }

# Make subplots and iterate over dictionary for hydrostatic imbalances
fig,axs = plt.subplots(nrows = 1, ncols = 4,figsize=(50,5),subplot_kw={'xticks':[],
                       'yticks':[]})
fig.suptitle('Deviation of hydrostatic equilibrium with channel widths 150, 250, 300 shift and 250', 
             fontsize = 30)

for ax, run in zip(axs.flat, hydrostatic_thicknesses.keys()):
    x = hydrostatic_thicknesses[run][0]
    y = hydrostatic_thicknesses[run][1]
    ht = hydrostatic_thicknesses[run][2]
    #thick_calc = hydrostatic_thicknesses[run][3]
    #thick_model = hydrostatic_thicknesses[run][4]
    im = ax.tripcolor(x,y,ht,shading='gouraud',vmin=-15,vmax=15,cmap = 'RdBu')   
    fig.colorbar(im, ax=ax)
    points = [x,y]
    points = np.asarray(points)
    points = points.transpose()
    
    
    # Delaunay triangulation for grid
    tri = Triangulation(x,y)
    delaunay = Delaunay(points)
    ax.triplot(x,y,delaunay.simplices,alpha=0.1)
    ax.title.set_text('width = {:s}, time = {:s}' .format(run.split('_')[2], run.split('_')[-1]))
   # ax.axis('equal')
   
   
   
   
   
   
   
   
   
   
   
   
   
neighbor = 1   
cs = 1070000                          
rcParams['figure.figsize'] = 50,20



x = ch[0]
y = ch[1]

x = np.ndarray.tolist(x)
y = np.ndarray.tolist(y)
orig_len = len(x)
x = x[-3:-1] + x + x[1:3]
y = y[-3:-1] + y + y[1:3]

t = np.arange(len(x))
ti = np.linspace(2, orig_len + 1, 10 * orig_len)

xi = interp1d(t, x, kind='cubic')(ti)
yi = interp1d(t, y, kind='cubic')(ti)

fig, ax = plt.subplots()
ax.plot(xi, yi)
ax.plot(x, y,'ro')
ax.margins(0.05)
plt.show()




plt.plot(x,y,'bo')
plt.show()

plt.plot(hullpoints[:,0],hullpoints[:,1],'ro')
plt.show()

   
   
   
   
   