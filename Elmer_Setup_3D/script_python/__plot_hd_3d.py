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
from __main__ import ModelRun
from scipy.spatial import Delaunay,ConvexHull
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
from math import trunc








class Plot_hydrostatic_deviation():
    
    """
    Class Plot_hydrostatic_deviation:
    Here description
    """
    

    
    
    
    def __init__(self,
                 t1 = None,
                 t2 = None,
                 t3 = None,
                 t4 = None,
                 amp = None,
                 width = None,
                 **kwargs): 
                             
                 
                
        """
        Plot_hydrostatic_deviation 
        
        Parameters
        ----------
        t1 : int
            timestep first plot
        t2 : int
            timestep second plot
        t3 : int
            timestep third plot  
        t4 : int
            timestep fourth plot
        amp : int
            amplitude of modelrun (see dictionary of modelrun-class
            to choose amplitude)
        width : int
            width of gauss function (valid for both directions)    
        """
        
        
        
        #======================================================================
        # Setup fonts
        #======================================================================
        

        font_annotation = {'color':'black',
                       'size': '17'
                }    
       
        font_title = {'color':'black',
               'size':'20',
               'weight':'bold'
               }

        font_axes = {'color':'black',
               'size':'14',
               'weight':'normal'
               }
        
        font_label = {'color':'black',
               'size':'15',
               'weight':'normal'
               }
        
        
        
        
        
        
        #======================================================================
        # Define variables for input
        #======================================================================
       
              
        # Choose default timesteps if timesteps are not given
        if  None is (t1,t2,t3,t4):
                t1 = 10
                t2 = 50
                t3 = 100
                t4 = 200
                
                
        self.t1 = t1
        self.t2 = t2
        self.t3 = t2
        self.t4 = t3
        self.amp = amp
        self.width = width
        
        # Crosssections
        cs1 = 1062500
        cs2 = 1070000
        cs3 = 1075500
        
        # Bounding box for zoom
        x1,x2,y1,y2 = cs1,cs2,-200,200
        
        # Boundaries for model domain
        ymin = -5000
        ymax = 5000
        GL = 1056000
        
        
        
        
        #======================================================================
        # Define figure and plot properties
        #======================================================================
        
        fig = plt.figure(figsize = (40,40))    
        
        gs_0 = gridspec.GridSpec(2,2, hspace = 0.2, wspace = 0.1, figure = fig)
        self.times = [t1, t2, t3, t4]
        
        # 2,2 plots with 3 subplots
        for i ,t in zip(range(4), self.times):
            gs00 = gridspec.GridSpecFromSubplotSpec(3,3, \
            height_ratios=[0.05,1,0.5], hspace = 0.25 ,subplot_spec=gs_0[i])
        
        # Run ModelRun-Class with default or self
            if None is (self.amp,self.width):
                mr = ModelRun(250,150,150,0,t)
            mr = ModelRun(self.amp,self.width,self.width,0,t)
            
            ht = mr.compute_hydrostatic_thickness()
            c1 = mr.compute_concavehull(cs1)
            c2 = mr.compute_concavehull(cs2)
            c3 = mr.compute_concavehull(cs3)
           
            
            x = ht[0]
            y = ht[1]
            ht_array =ht[2] 
            
            # Points for triangulation
            points = [x,y]
            points = np.asarray(points)
            points = points.transpose()
            
            
            
            #------------------------------------------------------------------
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
            original = c1[2]
            ax2.plot(lower,'b-', label = 'Modelled thickness')
            ax2.plot(upper,'b')
            ax2.plot(original,'r--', label = 'Hydrostatic thickness')
            ax2.legend(loc = 'lower right', bbox_to_anchor=(0.55,0.77))
            ax2.grid()
            ax2.set_xlabel('Width [m]', fontdict = font_label)
            ax2.set_ylabel('Height [m]', fontdict = font_label)
            ax2.set_ylim([-350,55])
            bottom, top = ax2.get_ylim()
            ax2.text(ymin,top,"A",fontdict=font_annotation,ha = 'right', va='bottom')
            ax2.text(ymax,top,"A'",fontdict=font_annotation,ha = 'left', va='bottom')
            fig.add_subplot(ax2)
           
#            a
            
            
            
            #------------------------------------------------------------------
            ax3 = plt.Subplot(fig,gs00[2,1])
            lower = c2[0]
            upper = c2[1]
            original = c2[2]
            ax3.plot(lower,'b-')
            ax3.plot (upper,'b')
            ax3.plot(original,'r--')
            
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
            original = c3[2]
            ax4.plot(lower,'b-')
            ax4.plot (upper,'b')
            ax4.plot(original,'r--')
            
            ax4.grid()
            ax4.set_ylim([-350,55])
            bottom, top = ax4.get_ylim()
            ax4.text(ymin, top,"C",fontdict=font_annotation,ha = 'right', va='bottom')
            ax4.text(ymax, top,"C'",fontdict=font_annotation,ha = 'left', va='bottom')
            fig.add_subplot(ax4)
            
        
                  
        fig.show()
        