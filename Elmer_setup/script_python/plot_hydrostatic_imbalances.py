#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 14:28:11 2019

@author: jloos
"""

import numpy as np
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
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
from scipy.interpolate import LinearNDInterpolator,Rbf, UnivariateSpline,CloughTocher2DInterpolator, griddata
import concave_hull
from collections import OrderedDict
from matplotlib.colors import LogNorm








# Get runs, put them into a dictionary (iterable)
hydrostatic_thicknesses= {
    '1_ht_250_25' : ModelRun(250,250,250,0,25).compute_hydrostatic_thickness(),
    '2_ht_250_50' : ModelRun(250,250,250,0,50).compute_hydrostatic_thickness(),
    '3_ht_250_100' : ModelRun(250,250,250,0,100).compute_hydrostatic_thickness(),
    '4_ht_250_150' : ModelRun(250,250,250,0,150).compute_hydrostatic_thickness(),
    
    '5_ht_500_25' : ModelRun(150,500,500,0,25).compute_hydrostatic_thickness(),
    '6_ht_500_50' : ModelRun(150,500,500,0,50).compute_hydrostatic_thickness(),
    '7_ht_500_100': ModelRun(150,500,500,0,100).compute_hydrostatic_thickness(),
    '8_ht_500_150': ModelRun(150,500,500,0,150).compute_hydrostatic_thickness()
    }

hydrostatic_thicknesses = OrderedDict(hydrostatic_thicknesses)




concave_hulls= {
    '1_ch_250_25' : ModelRun(250,250,250,0,25).compute_concavehull(1076000,2),
    '2_ch_250_50' : ModelRun(250,250,250,0,50).compute_concavehull(1076000,5),
    '3_ch_250_100' : ModelRun(250,250,250,0,100).compute_concavehull(1076000,5),
    '4_ch_250_150' : ModelRun(250,250,250,0,100).compute_concavehull(1076000,3),
    
    '5_ch_500_25' : ModelRun(150,500,500,0,25).compute_concavehull(1076000,3),
    '6_ch_500_50' : ModelRun(150,500,500,0,50).compute_concavehull(1076000,3),
    '7_ch_500_100': ModelRun(150,500,500,0,100).compute_concavehull(1076000,3),
    '8_ch_250_150' : ModelRun(250,250,250,0,100).compute_concavehull(1076000,3)
    }

concave_hulls = OrderedDict(concave_hulls)





# Define properties for subplots
font_title ={'color':'black',
       'size':'40',
       'weight':'bold'
       }

font_axes ={'color':'black',
       'size':'14',
       'weight':'italic'
       }

# Make subplots and iterate over dictionary for hydrostatic imbalances
fig,axs = plt.subplots(nrows = 2, ncols = 4,figsize=(40,15),subplot_kw={'xticks':[],
                       'yticks':[]})
fig.suptitle('Deviation of hydrostatic equilibrium with channel width 500 and 250', 
             fontsize = 30)

for ax, run in zip(axs.flat, hydrostatic_thicknesses.keys()):
    x = hydrostatic_thicknesses[run][0]
    y = hydrostatic_thicknesses[run][1]
    ht = hydrostatic_thicknesses[run][2]
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
    
    
    stress = ModelRun(250,250,250,0,50).cut_and_slice(1076000,'stress vector')[0] 
    linepoints = ModelRun(250,250,250,0,50).cut_and_slice(1076000,'stress vector')[1]    
    ax.plot(linepoints[:,0],linepoints[:,1],'-r',label = 'A')
    ax.text(1076000,-5000,'A',fontdict=font)
    ax.text(1076000,5000,"A'",fontdict=font)
    

        






# Make subplots and iterate over dictionary for concave hulls 
fig1,axs1 = plt.subplots(nrows = 2, ncols = 3,figsize=(30,15))
fig1.suptitle('Concave hulls of crosssections', 
             fontsize = 30)
    
for ax, run in zip(axs1.flat, concave_hulls.keys()):
    x = concave_hulls[run][:,0]
    z = concave_hulls[run][:,1]
   
    ax.plot(x,z,'b-')
    ax.grid()
    ax.title.set_text('width = {:s}, time = {:s}' .format(run.split('_')[2], run.split('_')[-1]))
    
   
x_list = linepoints[:,1]
y_list = linepoints[:,2]
z_list = stress[:,2]    
 
    
 
    
X, Y = np.meshgrid(x_list,y_list)

# Show the positions of the sample points, just to have some reference

#plt.scatter(x_list,y_list,400,facecolors='none')  
    
npoints = np.delete(linepoints, 0, 1)   
grid = griddata(npoints,z_list,(X,Y),method = 'nearest')  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
## Manually    
#ax = axs[0,0]
#ax.triplot(x,y,delaunay.simplices,alpha=0.2)
#im = ax.tripcolor(tri,ht_250_25[2],shading='gouraud',vmin=-5,vmax=15)
#fig.colorbar(im, ax=ax)
#ax.set_title("t=25")
#
#ax = axs[0,1]
#ax.triplot(x,y,delaunay.simplices,alpha=0.2)
#im = ax.tripcolor(tri,ht_500_25[2],shading='gouraud',vmin=-5,vmax=15)
#fig.colorbar(im, ax=ax)
#ax.set_title("t=25")
#
#
#ax = axs[1,0]
#ax.triplot(x,y,delaunay.simplices,alpha=0.2)
#im = ax.tripcolor(tri,ht_250_50[2],shading='gouraud',vmin=-5,vmax=15)
#fig.colorbar(im, ax=ax)
#ax.set_title("t=50")
#
#ax = axs[1,1]
#ax.triplot(x,y,delaunay.simplices,alpha=0.2)
#im = ax.tripcolor(tri,ht_500_50[2],shading='gouraud',vmin=-5,vmax=15)
#fig.colorbar(im, ax=ax)
#ax.set_title("t=50")
#
#
#ax = axs[2,0]
#ax.triplot(x,y,delaunay.simplices,alpha=0.2)
#im = ax.tripcolor(tri,ht_250_100[2],shading='gouraud',vmin=-5,vmax=15)
#fig.colorbar(im, ax=ax)
#ax.set_title("t=100")
#
#ax = axs[2,1]
#ax.triplot(x,y,delaunay.simplices,alpha=0.2)
#im = ax.tripcolor(tri,ht_500_100[2],shading='gouraud',vmin=-5,vmax=15)
#fig.colorbar(im, ax=ax)
#ax.set_title("t=100")
#plt.show()
