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
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar








#hydrostatic_thicknesses = {
#    '1_ht_150_10' : ModelRun(250,400,400,'double',t1).compute_hydrostatic_thickness(),
#    '2_ht_150_10' : ModelRun(250,400,400,'double',t2).compute_hydrostatic_thickness(),
#    '3_ht_150_10' : ModelRun(250,400,400,'double',t3).compute_hydrostatic_thickness(),
#    '4_ht_150_10' : ModelRun(250,400,400,'double',t4).compute_hydrostatic_thickness(),
#    
#                            }
#    
#hydrostatic_thicknesses = OrderedDict(hydrostatic_thicknesses)
#
#
#ch_1= {
#    '1.1_ch_150_5' : ModelRun(250,150,150,0,t1).compute_concavehull(cs1,neighbor),
#    '1.2_ch_150_5' : ModelRun(250,150,150,0,t1).compute_concavehull(cs2,neighbor),
#    '1.3_ch_150_5' : ModelRun(250,150,150,0,t1).compute_concavehull(cs3,neighbor),
#    
#    '2.1_ch_150_5' : ModelRun(250,150,150,0,t2).compute_concavehull(cs1,neighbor),
#    '2.2_ch_150_5' : ModelRun(250,150,150,0,t2).compute_concavehull(cs2,neighbor),
#    '2.3_ch_150_5' : ModelRun(250,150,150,0,t2).compute_concavehull(cs3,neighbor),
#    
#    '3.1_ch_150_5' : ModelRun(250,150,150,0,t3).compute_concavehull(cs1,neighbor),
#    '3.2_ch_150_5' : ModelRun(250,150,150,0,t3).compute_concavehull(cs2,neighbor),
#    '3.3_ch_150_5' : ModelRun(250,150,150,0,t3).compute_concavehull(cs3,neighbor),
#    
#    '4.1_ch_150_5' : ModelRun(250,150,150,0,t4).compute_concavehull(cs1,neighbor),
#    '4.2_ch_150_5' : ModelRun(250,150,150,0,t4).compute_concavehull(cs2,neighbor),
#    '4.3_ch_150_5' : ModelRun(250,150,150,0,t4).compute_concavehull(cs3,neighbor),
#    
#    
#    }

#for t in times:
#    mr = ModelRun(250,150,150,0,t)
#    
#    ht = mr.compute_hydrostatic_thickness()
#    c1 = mr.compute_concavehull(cs1,neighbor)
#    c2 = mr.compute_concavehull(cs2,neighbor)
#    c3 = mr.compute_concavehull(cs3,neighbor)
#    
#    
#ch_1 = OrderedDict(ch_1)
#
#
#
#ax1 = plt.subplot2grid((3,3), (0,0), colspan=3)
#ax3 = plt.subplot2grid((3,3), (1, 0))
#ax4 = plt.subplot2grid((3,3), (1, 1))
#ax5 = plt.subplot2grid((3,3), (1, 2))


 # Make subplots and iterate over dictionary for hydrostatic imbalances
#fig,axs = plt.subplots(nrows = 4, ncols = 4,figsize=(50,20))
#fig.suptitle('Deviation of hydrostatic equilibrium with channel widths 150, 250, 300 shift and 250', 
 #            fontsize = 30)
 
 
#ax1 = plt.subplot2grid((3,3), (0,0), colspan=3)
#ax3 = plt.subplot2grid((3,3), (1, 0))
#ax4 = plt.subplot2grid((3,3), (1, 1))
#ax5 = plt.subplot2grid((3,3), (1, 2)) 



for ax, run in zip(ax1, hydrostatic_thicknesses.keys()):
    x = hydrostatic_thicknesses[run][0]
    y = hydrostatic_thicknesses[run][1]
    ht = hydrostatic_thicknesses[run][2]
    im = ax1.tripcolor(x,y,ht,shading='gouraud',vmin=-15,vmax=15,cmap = 'RdBu')
    cbar = fig.colorbar(im, ax1=ax1, extend = 'both')
    cbar.set_label('Deviation of hydrostatic equilibrium [m]')
    axins = zoomed_inset_axes(ax1,2,loc='upper right')
    axins.tripcolor(x,y,ht,shading='gouraud',vmin=-15,vmax=15,cmap = 'RdBu')
    axins.set_xlim(x1,x2)
    axins.set_ylim(y1,y2)
    plt.yticks(visible=False)
    plt.xticks(visible=False)
    mark_inset(ax1,axins,loc1=2,loc2=4,fc="none", ec="0.5")


#========================================
#========================================
font_annotation = {'color':'black',
                   'size': '10'
                   }    
    
    
    
x1,x2,y1,y2 = 1070000,1075000,-1000,1000
t1 = 10
t2 = 50
t3 = 100
t4 = 200
cs1 = 1062500
cs2 = 1070000
cs3 = 1075500
ymin = -5000
ymax = 5000

neighbor = 4



fig = plt.figure(figsize = (60,10))    

gs_0 = gridspec.GridSpec(1,4, hspace = 0.1, wspace = 0.1, figure = fig)
times = [t1, t2, t3, t4]

# 1,4 plots with 3 subplots
for i ,t in zip(range(4), times):
    gs00 = gridspec.GridSpecFromSubplotSpec(3,3, height_ratios=[0.05,1,0.5], hspace = 0.15 ,subplot_spec=gs_0[i])
    
    mr = ModelRun(250,150,150,0,t)
    
    ht = mr.compute_hydrostatic_thickness()
    c1 = mr.compute_concavehull(cs1,neighbor)
    c2 = mr.compute_concavehull(cs2,neighbor)
    c3 = mr.compute_concavehull(cs3,neighbor)
   
    x = ht[0]
    y = ht[1]
    ht_array =ht[2] 
    
    # Points for triangulation
    points = [x,y]
    points = np.asarray(points)
    points = points.transpose()
    
    
    
    #-------------------------------
    ax1 = plt.Subplot(fig,gs00[1,:])
    im = ax1.tripcolor(x,y,ht_array,shading='gouraud',vmin=-15,vmax=15,cmap = 'RdBu')
    
    # Colorbar definition (extra axis)
    cbar = plt.subplot(gs00[0,:])
    cbar = Colorbar(ax=cbar,mappable = im, extend = 'both', orientation = 'horizontal', ticklocation = 'top')
    cbar.set_label('Deviation of hydrostatic equilibrium [m]',labelpad=10)
    
 
    # Crosssections
    line_cs1 = ax1.axvline(x=cs1, ymin=ymin, ymax=ymax, color='r')
    line_cs2 = ax1.axvline(x=cs2, ymin=ymin, ymax=ymax, color='r')
    line_cs3 = ax1.axvline(x=cs3, ymin=ymin, ymax=ymax, color='r')
    
    ax1.text(cs1,-5000,"A",fontdict=font_annotation)
    ax1.text(cs1,5000,"A'",fontdict=font_annotation)
    ax1.text(cs2,-5000,"B",fontdict=font_annotation)
    ax1.text(cs2,5000,"B'",fontdict=font_annotation)
    ax1.text(cs3,-5000,"C",fontdict=font_annotation)
    ax1.text(cs3,5000,"C'",fontdict=font_annotation)
    
    
    # Delaunay triangulation for grid
    tri = Triangulation(x,y)
    delaunay = Delaunay(points)
    ax1.triplot(x,y,delaunay.simplices,alpha=0.05)
    
    fig.add_subplot(ax1)
    
    # Zoomed in rectangle for better visualisation of channel
    axins = zoomed_inset_axes(ax1,2,loc='upper right')
    axins.tripcolor(x,y,ht_array,shading='gouraud',vmin=-15,vmax=15,cmap = 'RdBu')
    axins.set_xlim(x1,x2)
    axins.set_ylim(y1,y2)
    plt.yticks(visible=False)
    plt.xticks(visible=False)
    mark_inset(ax1,axins,loc1=2,loc2=4,fc="none", ec="0.5")
        
   
    
    #-------------------------------
    ax2 = plt.Subplot(fig,gs00[2,0])
    lower = c1[0]
    upper = c1[1]       
    ax2.plot(lower,'b-')
    ax2.plot (upper,'b')
    ax2.grid()
    fig.add_subplot(ax2)
   
    #-------------------------------
    ax3 = plt.Subplot(fig,gs00[2,1])
    lower = c2[0]
    upper = c2[1]       
    ax3.plot(lower,'b-')
    ax3.plot (upper,'b')
    ax3.grid()
    fig.add_subplot(ax3)
    
    #-------------------------------
    ax4 = plt.Subplot(fig,gs00[2,2])
    lower = c3[0]
    upper = c3[1]       
    ax4.plot(lower,'b-')
    ax4.plot (upper,'b')
    ax4.grid()
    fig.add_subplot(ax4)
    

          
fig.show()







