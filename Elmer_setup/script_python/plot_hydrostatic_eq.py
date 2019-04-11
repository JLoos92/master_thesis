#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 15:29:02 2019

@author: jloos
"""

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
from pylab import rcParams
import numpy as np
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
import turbulucid
from modelrun import ModelRun





class Plot_hydrostatic_deviation():
    
    """
    Class Plot:
    Here description
    """
    
# Check if volume is mounted, if not mount volume
    
    
    
    def __init__(self,
                 **kwargs): 
                             
                 
                
        """
        Plot_hydrostatic_deviation 
        
        Parameters
        ----------
       
        """
        
        # Setup fonts
        

        font_annotation = {'color':'black',
                       'size': '20'
                       }    
       
        font_title = {'color':'black',
               'size':'40',
               'weight':'bold'
               }

        font_axes = {'color':'black',
               'size':'14',
               'weight':'italic'
               }
        
        # define variables
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
        self.times = [t1, t2, t3, t4]
        
        # 1,4 plots with 3 subplots
        for i ,t in zip(range(4), self.times):
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
            cbar = Colorbar(ax=cbar, mappable = im, extend = 'both', orientation = 'horizontal', ticklocation = 'top')
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
