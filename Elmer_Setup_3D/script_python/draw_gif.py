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
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar    

def plot(t):

     #======================================================================
     # Setup fonts
     #======================================================================
        

    font_annotation = {'color':'black',
                   'size': '16'
            }    
   
    font_title = {'color':'black',
           'size':'22',
           'weight':'bold'
           }

    font_axes = {'color':'black',
           'size':'14',
           'weight':'italic'
           }
    
    font_label = {'color':'black',
           'size':'16',
           'weight':'normal'
           }
    
        
    # Crosssections
    cs1 = 1062500
    cs2 = 1070000
    cs3 = 1075500   
        
        
        
    # Bounding box for zoom
    x1,x2,y1,y2 = cs2,cs3,-1000,1000
    
    # Boundaries for model domain
    ymin = -5000
    ymax = 5000
    GL = 1056000
    
    # Probably not necessary, check in modelrun class
    neighbor = 4 
    
    fig = plt.figure(figsize = (30,20))    
    
    gs_0 = gridspec.GridSpec(1,1, hspace = 0.1, wspace = 0.1, figure = fig)
    
    

    gs00 = gridspec.GridSpecFromSubplotSpec(3,3, height_ratios=[0.05,1,0.5], hspace = 0.25 ,subplot_spec=gs_0[0])
    
      
    mr = ModelRun(250,100,100,0,t)
            
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
    cbar = Colorbar(ax=cbar, mappable = im, extend = 'both', \
    orientation = 'horizontal', ticklocation = 'top') 
    cbar.set_label('Deviation of hydrostatic equilibrium [m] \n' \
    + 't = ' + str(t) + ', a = ' + str(self.amp) + ', width = ' \
    + str(self.width) ,labelpad=10, fontdict = font_title)
    
    ax1.set_xlabel('Distance from grounding line [m]',fontdict = font_label)
    ax1.set_ylabel('Along flow direction [m]',fontdict = font_label)
    start,end = ax1.get_xlim()
    ax1.set_xticklabels(list(map(int,(np.arange(start,end, 2500)-GL))))
   
    
 
    # Crosssections
    line_cs1 = ax1.vlines(x=cs1, ymin=ymin, ymax=ymax, color='r')
    line_cs2 = ax1.vlines(x=cs2, ymin=ymin, ymax=ymax, color='r')
    line_cs3 = ax1.vlines(x=cs3, ymin=ymin, ymax=ymax, color='r')
    
    #s,e = yaxis = ax1.get_ylim()
    
    ax1.text(cs1,-5000,"A",fontdict=font_annotation,ha = 'center', va='top')
    ax1.text(cs1,5000,"A'",fontdict=font_annotation,ha = 'center', va='bottom')
    ax1.text(cs2,-5000,"B",fontdict=font_annotation,ha = 'center', va='top')
    ax1.text(cs2,5000,"B'",fontdict=font_annotation,ha = 'center', va='bottom')
    ax1.text(cs3,-5000,"C",fontdict=font_annotation,ha = 'center', va='top')
    ax1.text(cs3,5000,"C'",fontdict=font_annotation,ha = 'center', va='bottom')
    
    
    # Delaunay triangulation for grid
    delaunay = Delaunay(points)
    ax1.triplot(x,y,delaunay.simplices,alpha=0.05)
    
    fig.add_subplot(ax1)
    
    # Zoomed in rectangle for better visualisation of channel
    axins = zoomed_inset_axes(ax1,1.5,loc='upper right',borderpad=5)
    axins.tripcolor(x,y,ht_array,shading='gouraud',vmin=-15,vmax=15,cmap = 'RdBu')
    axins.set_xlim(x1,x2)
    axins.set_ylim(y1,y2)
    plt.yticks(visible=False)
    plt.xticks(visible=False)
    mark_inset(ax1,axins,loc1=2,loc2=4,fc="none", ec="0.5")
        
   
    
    #------------------------------------------------------------------
    ax2 = plt.Subplot(fig,gs00[2,0])
    lower = c1[0]
    upper = c1[1]       
    ax2.plot(lower,'b-')
    ax2.plot (upper,'b')
    ax2.grid()
    ax2.set_xlabel('Width [m]', fontdict = font_label)
    ax2.set_ylabel('Height [m]', fontdict = font_label)
    ax2.set_ylim([-350,55])
    bottom, top = ax2.get_ylim()
    ax2.text(ymin,top,"A",fontdict=font_annotation,ha = 'right', va='bottom')
    ax2.text(ymax,top,"A'",fontdict=font_annotation,ha = 'left', va='bottom')
    fig.add_subplot(ax2)
   
    #------------------------------------------------------------------
    ax3 = plt.Subplot(fig,gs00[2,1])
    lower = c2[0]
    upper = c2[1]       
    ax3.plot(lower,'b-')
    ax3.plot (upper,'b')
    ax3.grid()
    ax3.set_ylim([-350,55])
    bottom, top = ax3.get_ylim()
    ax3.text(ymin,top,"B",fontdict=font_annotation,ha = 'right', va='bottom')
    ax3.text(ymax,top,"B'",fontdict=font_annotation,ha = 'left', va='bottom')
    fig.add_subplot(ax3)
    
    #------------------------------------------------------------------
    ax4 = plt.Subplot(fig,gs00[2,2])
    lower = c3[0]
    upper = c3[1]       
    ax4.plot(lower,'b-')
    ax4.plot (upper,'b')
    ax4.grid()
    ax4.set_ylim([-350,55])
    bottom, top = ax4.get_ylim()
    ax4.text(ymin, top,"C",fontdict=font_annotation,ha = 'right', va='bottom')
    ax4.text(ymax, top,"C'",fontdict=font_annotation,ha = 'left', va='bottom')
    fig.add_subplot(ax4)
    
    


    fig.canvas.draw()       # draw the canvas, cache the renderer
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return image


kwargs_write = {'fps':5.0, 'quantizer':'nq'}
imageio.mimsave('./animation.gif', [plot(t) for t in range(100)], fps=5)

