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
from main import ModelRun
from scipy.spatial import Delaunay,ConvexHull
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
from math import trunc
from __plot_params import params_3d,params_horizontal

import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'






class Plot_hydrostatic_deviation():
    
    """
    Class Plot_hydrostatic_deviation:
        
    Plots the hydrostatic deviation at four timesteps which are defined by 
    t1 - t4. The amplitude of the bedrock bump in the 3D-case is not changed
    and will remain at 250m, whereas the hydrostatic deviation is a function of
    the channel width, the user can change the width of the bump for this spec-
    ific class. The plot shows a figure which contains one overview plot with
    three crosssections which can be altered. For each crosssection one subplot 
    is generated which shows the modelled hydrostatic thickness and calculated
    hydrostatic thickness at the vertical crosssection.
    """
    

    
    
    
    def __init__(self,
                 t1 = None,
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
        # Define variables for input
        #======================================================================
       
              
        # Choose default timesteps if timesteps are not given
        if  None is (t1):
                t1 = 10
                
                
                
        self.t1 = t1
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
        ymin_2 = -2000
        ymax_2 = 2000
        GL = 1056000
        
        # Colors
        orange = '#D55E00'
        blue = '#396AB1'    
        red = '#CC2529'
        wheat = '#948B3D'
        
        
        
        #======================================================================
        # Define figure and plot properties
        #======================================================================
        
        fig = plt.figure()    
        
        gs_0 = gridspec.GridSpec(1,1, wspace = 0, figure = fig)
        self.times = [t1]
        plt.rcParams.update(params_horizontal)
        fig.tight_layout()
        # 2,2 plots with 3 subplots
        for i ,t in zip(range(4), self.times):
            gs00 = gridspec.GridSpecFromSubplotSpec(3,3, \
            height_ratios=[0.05,0.8,0.5],subplot_spec=gs_0[i])
        
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
            cbar.set_label('Hydrostatic dev. [m]',labelpad=0)
            cbar.ax.tick_params(direction='in',length=1,width=1,pad=0)
            
            
            ax1.set_xlabel('Distance from grounding line [m]')
            ax1.set_ylabel('Across flow direction [m]')
            start,end = ax1.get_xlim()
            #ax1.set_xticklabels(list(map(int,(np.arange(start,end, 2500)-GL))))
            ax1.xaxis.tick_top()
            #ax1.xaxis.set_major_locator(ticker.MaxNLocator(3))
            ax1.xaxis.set_tick_params(pad=0)
            ax1.set_xticks([cs1,cs2,cs3])
            ax1.set_xticklabels([cs1-GL,cs2-GL,cs3-GL])
            
            ax1.spines['top'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['left'].set_visible(False)
            #ax1.tick_params(direction='in',length=3,width=1)
         
            # Crosssections
            line_cs1 = ax1.vlines(x=cs1, ymin=ymin, ymax=ymax, color='r')
            line_cs2 = ax1.vlines(x=cs2, ymin=ymin, ymax=ymax, color='r')
            line_cs3 = ax1.vlines(x=cs3, ymin=ymin, ymax=ymax, color='r')
            
            #s,e = yaxis = ax1.get_ylim()
            
            ax1.text(cs1,-5000,"A",ha = 'center', va='bottom')
            ax1.text(cs1,5000,"A'",ha = 'center', va='top')
            ax1.text(cs2,-5000,"B",ha = 'center', va='bottom')
            ax1.text(cs2,5000,"B'",ha = 'center', va='top')
            ax1.text(cs3,-5000,"C",ha = 'center', va='bottom')
            ax1.text(cs3,5000,"C'",ha = 'center', va='top')
            
            plt.setp(ax1.get_xticklabels(),fontweight = 'bold')
            plt.setp(ax1.get_yticklabels(),fontweight = 'bold')
            
            # Delaunay triangulation for grid
            delaunay = Delaunay(points)
            ax1.triplot(x,y,delaunay.simplices,alpha=0.05)
            ax1.set_yticks([-4000,-2000,0,2000,4000])
            ax1.set_ylim(-5000,5000)
            ax1.set_xlim(1060000,1079000)
            #ax1.autoscale_view('tight')
            plt.subplots_adjust(hspace=0.42)
            
            
            # place text box in upper left in axes coords      
            props = dict(boxstyle='round', facecolor='wheat')
            ax1.text(0.2, 0.92,'(b) t = ' + str(t1) + 'a\n @ cw = ' + str(self.width*2) + 'm', transform=ax1.transAxes, 
                     verticalalignment='top', bbox=props,weight='bold')
            
            
            fig.add_subplot(ax1)
            
         
            #------------------------------------------------------------------
            ax2 = plt.Subplot(fig,gs00[2,0])
            upper = c1[0]
            lower_2 = c1[1]
            original_2 = c1[2]
            ax2.plot(lower_2,'k-')
            #ax2.plot(upper,'b')
            ax2.plot(original_2,linestyle = '--', color = blue)
#            legend_ax2 = ax2.legend(['Modelled thickness','Hydrostatic thickness'],loc="upper right", prop=dict(weight='bold'))
#            frame_ax2 = legend_ax2.get_frame()
#            frame_ax2.set_facecolor('0.7')
#            frame_ax2.set_edgecolor('0.7')
            
            ax2.set_ylabel('Shelf elevation [m]')
            ax2.set_xlim(-2000,2000)
            ax2.set_ylim([-350,-150])
            ax2.set_yticklabels([-300,-250,-200])
            bottom, top = ax2.get_ylim()
            ax2.text(ymin_2,top,"A",ha = 'right', va='bottom')
            ax2.text(ymax_2,top,"A'",ha = 'left', va='bottom')
            ax2.set_xticks([-1500,0,1500])
            fig.add_subplot(ax2)

            
            #------------------------------------------------------------------
            ax3 = plt.Subplot(fig,gs00[2,1])
            upper = c2[0]
            lower_3 = c2[1]
            original_3 = c2[2]
            ax3.plot(lower_3,'k-')
            #ax3.plot (upper,'b')
            ax3.plot(original_3,linestyle = '--', color = blue)
            
            ax3.set_xlim(-2000,2000)
            ax3.set_ylim([-350,-150])
            ax3.set_yticks([-300,-250,-200])
            bottom, top = ax3.get_ylim()
            ax3.set_xlabel('Width [m]')
            ax3.text(ymin_2,top,"B",ha = 'right', va='bottom')
            ax3.text(ymax_2,top,"B'",ha = 'left', va='bottom')
            ax3.axes.get_yaxis().set_visible(False)
            ax3.set_xticks([-1500,0,1500])
            fig.add_subplot(ax3)
            
            #------------------------------------------------------------------
            ax4 = plt.Subplot(fig,gs00[2,2])
            upper = c3[0]
            lower_4 = c3[1]
            original_4 = c3[2]
            ax4.plot(lower_4,'k-')
            #ax4.plot (upper,'b')
            ax4.plot(original_4,linestyle = '--', color = blue)
            
            ax4.set_xlim(-2000,2000)
            ax4.set_ylim([-350,-150])
            ax3.set_yticks([-300,-250,-200])
            bottom, top = ax4.get_ylim()
            ax4.text(ymin_2, top,"C",ha = 'right', va='bottom')
            ax4.text(ymax_2, top,"C", ha = 'left', va='bottom')
            ax4.axes.get_yaxis().set_visible(False)
            ax4.set_xticks([-1500,0,1500])
            
            fig.add_subplot(ax4)
            #------------------------------------------------------------------
            
            ax4_twin = ax4.twinx()
            ax4_twin.plot(original_4-lower_4,color=red,linestyle=':')
            ax4_twin.set_ylim(-15,15)
            ax4_twin.set_yticks([-10,-5,0,5,10])
            ax4_twin.set_ylabel('HD [m]',color=red)
            
            ax3_twin = ax3.twinx()
            ax3_twin.plot(original_3-lower_3,color=red,linestyle=':')
            ax3_twin.set_yticks([-10,-5,0,5,10])
            ax3_twin.set_yticklabels([])
            ax3_twin.set_ylim(-15,15)
            
            ax2_twin = ax2.twinx()
            ax2_twin.plot(original_2-lower_2,color=red,linestyle=':')
            ax2_twin.set_yticks([-10,-5,0,5,10])
            ax2_twin.set_yticklabels([])
            ax2_twin.set_ylim(-15,15)
             
            plt.setp(ax4_twin.get_yticklabels(),color=red)
            
            
        path = str('plots/')
        fname= str('3d_' + str(self.amp) + str(self.width)+ '_' + str(t1) + '.png')
        fname_pdf = str('3d_' + str(self.amp) + str(self.width)+ '_' + str(t1) + '.pdf')
        
        fig.savefig(path + fname, format = 'png', dpi=1000)
        fig.savefig(path + fname_pdf, format = 'pdf', dpi=1000)
                  
        plt.show()
        