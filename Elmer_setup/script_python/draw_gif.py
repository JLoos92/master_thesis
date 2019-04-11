#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 10:33:04 2019

@author: jloos
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import gridspec
import numpy as np
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
import turbulucid
from modelrun import ModelRun
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
import imageio
        

def plot(t):

    font_annotation = {'color':'black',
                       'size': '10'
                       }    
        
        
        
    x1,x2,y1,y2 = 1070000,1075000,-1000,1000
    
    cs1 = 1062500
    cs2 = 1070000
    cs3 = 1075500
    ymin = -5000
    ymax = 5000
    
    neighbor = 4
    
    

    fig = plt.figure(figsize = (30,15))    
    
    
    
    gs_0 = gridspec.GridSpec(1,1, hspace = 0.1, wspace = 0.1, figure = fig)
    
    

    gs00 = gridspec.GridSpecFromSubplotSpec(3,3, height_ratios=[0.05,1,0.5], hspace = 0.15 ,subplot_spec=gs_0[0])
    
      
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
    im = ax1.tripcolor(x,y,ht_array,shading='gouraud',vmin=-15,vmax=15,cmap=plt.cm.get_cmap('RdBu',101))
    
    # Colorbar definition (extra axis)
    cbar = plt.subplot(gs00[0,:])
    cbar = Colorbar(ax=cbar,mappable = im, extend = 'both', orientation = 'horizontal', ticklocation = 'top')
    cbar.set_label('Deviation of hydrostatic equilibrium [m]',labelpad=10)
    cbar.set_clim([-15,15])
     
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
    axins.tripcolor(x,y,ht_array,shading='gouraud',vmin=-15,vmax=15,cmap=plt.cm.get_cmap('RdBu',101))
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
    ax2.set_ylim([-350,50])
    fig.add_subplot(ax2)
       
    #-------------------------------
    ax3 = plt.Subplot(fig,gs00[2,1])
    lower = c2[0]
    upper = c2[1]       
    ax3.plot(lower,'b-')
    ax3.plot (upper,'b')
    ax3.grid()
    ax3.set_ylim([-350,50])
    fig.add_subplot(ax3)
    
    #-------------------------------
    ax4 = plt.Subplot(fig,gs00[2,2])
    lower = c3[0]
    upper = c3[1]       
    ax4.plot(lower,'b-')
    ax4.plot (upper,'b')
    ax4.grid()
    ax4.set_ylim([-350,50])
    fig.add_subplot(ax4)     


    fig.canvas.draw()       # draw the canvas, cache the renderer
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return image


kwargs_write = {'fps':4.0, 'quantizer':'nq'}
imageio.mimsave('./powers.gif', [plot(t) for t in range(200)], fps=4)

