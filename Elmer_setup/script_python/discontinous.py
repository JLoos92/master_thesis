#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 16:10:11 2019

@author: jloos
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
from scipy.interpolate import griddata,interp2d
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


t1 = 10
t2 = 50
t3 = 100
t4 = 200

neighbor = 4



# Get runs, put them into a dictionary (iterable)
hydrostatic_thicknesses = {
    '1.1_ht_150_10' : ModelRun(250,150,150,0,t1).compute_hydrostatic_thickness(),
    '1.2_ht_150_25' : ModelRun(250,150,150,0,t2).compute_hydrostatic_thickness(),
    '1.3_ht_150_50' : ModelRun(250,150,150,0,t3).compute_hydrostatic_thickness(),
    '1.4_ht_150_100' : ModelRun(250,150,150,0,t4).compute_hydrostatic_thickness(),
    
    '2.1_ht_150shift_10' : ModelRun(250,250,250,0,t1).compute_hydrostatic_thickness(),
    '2.2_ht_150shift_25' : ModelRun(250,250,250,0,t2).compute_hydrostatic_thickness(),
    '2.3_ht_150shift_50': ModelRun(250,250,250,0,t3).compute_hydrostatic_thickness(),
    '2.4_ht_150shift_100': ModelRun(250,250,250,0,t4).compute_hydrostatic_thickness(),
    
    '3.1_ht_300_10' : ModelRun(250,400,400,'double',t1).compute_hydrostatic_thickness(),
    '3.2_ht_300_25' : ModelRun(250,400,400,'double',t2).compute_hydrostatic_thickness(),
    '3.3_ht_300_50': ModelRun(250,400,400,'double',t3).compute_hydrostatic_thickness(),
    '3.4_ht_300_100': ModelRun(250,400,400,'double',t4).compute_hydrostatic_thickness(),
    
    '4.1_ht_400_10' : ModelRun(250,300,300,'double',t1).compute_hydrostatic_thickness(),
    '4.2_ht_400_25' : ModelRun(250,300,300,'double',t2).compute_hydrostatic_thickness(),
    '4.3_ht_400_50': ModelRun(250,300,300,'double',t3).compute_hydrostatic_thickness(),
    '4.4_ht_400_100': ModelRun(250,300,300,'double',t4).compute_hydrostatic_thickness()
    }

hydrostatic_thicknesses = OrderedDict(hydrostatic_thicknesses)

cs = 1060000

concave_hulls= {
    '1.1_ch_150_5' : ModelRun(250,150,150,0,t1).compute_concavehull(cs,neighbor),
    '1.2_ch_150_10' : ModelRun(250,150,150,0,t2).compute_concavehull(cs,neighbor),
    '1.3_ch_150_25' : ModelRun(250,150,150,0,t3).compute_concavehull(cs,neighbor),
    '1.4_ht_150_62' : ModelRun(250,150,150,0,t4).compute_concavehull(cs,neighbor),
    
    '2.1_ch_250_5' : ModelRun(250,250,250,0,t1).compute_concavehull(cs,neighbor),
    '2.2_ch_250_10' : ModelRun(250,250,250,0,t2).compute_concavehull(cs,neighbor),
    '2.3_ch_250_25': ModelRun(250,250,250,0,t3).compute_concavehull(cs,neighbor),
    '2.4_ch_250_62': ModelRun(250,250,250,0,t4).compute_concavehull(cs,neighbor),
    
    '3.1_ch_300_5' : ModelRun(250,400,400,'double',t1).compute_concavehull(cs,neighbor),
    '3.2_ch_300_10' : ModelRun(250,400,400,'double',t2).compute_concavehull(cs,neighbor),
    '3.3_ch_300_25': ModelRun(250,400,400,'double',t3).compute_concavehull(cs,neighbor),
    '3.4_ch_300_62': ModelRun(250,400,400,'double',t4).compute_concavehull(cs,neighbor),
    
    '4.1_ch_300_10' : ModelRun(250,300,300,'double',t1).compute_concavehull(cs,neighbor),
    '4.2_ch_300_25' : ModelRun(250,300,300,'double',t2).compute_concavehull(cs,neighbor),
    '4.3_ch_300_50': ModelRun(250,300,300,'double',t3).compute_concavehull(cs,neighbor),
    '4.4_ch_300_100': ModelRun(250,300,300,'double',t4).compute_concavehull(cs,neighbor),
    
    
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

font_annotation = {'color':'black',
                   'size': '10'
                   }
    
  

 # Make subplots and iterate over dictionary for hydrostatic imbalances
fig,axs = plt.subplots(nrows = 4, ncols = 4,figsize=(50,20))
fig.suptitle('Deviation of hydrostatic equilibrium with channel widths 150, 250, 300 shift and 250', 
             fontsize = 30)

for ax, run in zip(axs.flat, hydrostatic_thicknesses.keys()):
    x = hydrostatic_thicknesses[run][0]
    y = hydrostatic_thicknesses[run][1]
    ht = hydrostatic_thicknesses[run][2]
    im = ax.tripcolor(x,y,ht,shading='gouraud',vmin=-15,vmax=15,cmap = 'RdBu')
    cbar = fig.colorbar(im, ax=ax, extend = 'both')
    cbar.set_label('Deviation of hydrostatic equilibrium [m]')
    axins = zoomed_inset_axes(ax,2,loc='upper right')
    axins.tripcolor(x,y,ht,shading='gouraud',vmin=-15,vmax=15,cmap = 'RdBu')
    axins.set_xlim(x1,x2)
    axins.set_ylim(y1,y2)
    plt.yticks(visible=False)
    plt.xticks(visible=False)
    mark_inset(ax,axins,loc1=2,loc2=4,fc="none", ec="0.5")
    
    
    points = [x,y]
    points = np.asarray(points)
    points = points.transpose()
    
    
    # Delaunay triangulation for grid
    tri = Triangulation(x,y)
    delaunay = Delaunay(points)
    ax.triplot(x,y,delaunay.simplices,alpha=0.05)
    ax.title.set_text('width = {:s}, time = {:s}' .format(run.split('_')[2], run.split('_')[-1]))


  
#   # stress = ModelRun(250,150,150,0,50).cut_and_slice(1067000,'stress vector')[0] 
#    linepoints = ModelRun(250,150,150,0,10).cut_and_slice(1067000,'stress vector')[1]    
#    ax.plot(linepoints[:,0],linepoints[:,1],'-r',label = "A")
#    linepoints70 = ModelRun(250,150,150,0,10).cut_and_slice(1070000,'stress vector')[1]    
#    ax.plot(linepoints70[:,0],linepoints70[:,1],'-r',label = "B")
#    linepoints75 = ModelRun(250,150,150,0,10).cut_and_slice(1075000,'stress vector')[1]    
#    ax.plot(linepoints75[:,0],linepoints75[:,1],'-r',label = "C")
#    ax.text(1067000,-5000,"A",fontdict=font_annotation)
#    ax.text(1067000,5000,"A'",fontdict=font_annotation)
#    ax.text(1070000,-5000,"B",fontdict=font_annotation)
#    ax.text(1070000,5000,"B'",fontdict=font_annotation)
#    ax.text(1075000,-5000,"C",fontdict=font_annotation)
#    ax.text(1075000,5000,"C'",fontdict=font_annotation)
    

  
# Make subplots and iterate over dictionary for concave hulls 
fig1,axs1 = plt.subplots(nrows = 4, ncols = 4,figsize=(50,20))
fig1.suptitle('Concave hulls of crosssections at C', 
             fontsize = 30)
    
for ax, run in zip(axs1.flat, concave_hulls.keys()):
    lower = concave_hulls[run][1]
    upper = concave_hulls[run][2]
   
    ax.plot(lower,'b-')
    ax.plot (upper,'b')
    ax.grid()
    ax.title.set_text('width = {:s}, time = {:s}' .format(run.split('_')[2], run.split('_')[-1]))

