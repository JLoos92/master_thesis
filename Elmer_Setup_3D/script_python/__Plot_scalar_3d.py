#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 15:29:02 2019

@author: jloos
"""

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset, inset_axes
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from pylab import rcParams
import numpy as np
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
from main import ModelRun
from scipy.spatial import Delaunay,ConvexHull
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
from math import trunc
from __plot_params import params_3d,params_horizontal
import matplotlib.tri as tri


import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'






class Plot_scalar_3d():
    
    """
    Class: Plot_scalar_3d()
    -----------------------------------
        
    Plots the hydrostatic deviation for one timestep. 
    The amplitude of the bedrock bump in the 3D-case is not changed
    and will remain at 250m, whereas the hydrostatic deviation is a function of
    the channel width. The user can change the width of the bump for this spec-
    ific class. The plot shows a figure which contains a footprint plot with
    three crosssections which can be altered. This footprint plot is a 
    cropped (better visualisation) ice-shelf from the model output. For each 
    crosssection one subplot is generated which shows the modelled hydrostatic 
    thickness and calculated hydrostatic thickness at the vertical crosssection.
    Furthermore, the horizontal shear stress is depicted.
    """
    

    
    
    
    def __init__(self,
                 t1 = None,
                 amp = None,
                 width = None): 
                             
                 
                
        """
        Method : init 
        
        Parameters
        ----------
        t1 : int
            timestep first plot
        amp : int
            amplitude of modelrun (see dictionary of modelrun-class
            to choose amplitude)
        width : int
            width of gauss function (valid for both directions x,y)    
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
        ymin = -2000
        ymax = 2000
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
        
        # define gridspec
        gs_0 = gridspec.GridSpec(1,1, wspace = 0, figure = fig)
        self.times = [t1]
        plt.rcParams.update(params_horizontal)
        fig.tight_layout()
        
        
        # 2,2 plots with 3 subplots
        for i ,t in zip(range(4), self.times):
            gs00 = gridspec.GridSpecFromSubplotSpec(2,3, \
            height_ratios=[0.8,0.5],subplot_spec=gs_0[i],hspace=0.7)
        
            # Run ModelRun-Class with default or self
            if None is (self.amp,self.width):
                mr = ModelRun(250,150,150,0,t)
            mr = ModelRun(self.amp,self.width,self.width,'new',t)
            
            ht = mr.compute_hydrostatic_thickness()
            c1 = mr.compute_concavehull(cs1)
            c2 = mr.compute_concavehull(cs2)
            c3 = mr.compute_concavehull(cs3)

            scalar = mr.get_top_scalar('velocity')
            
            scalar_x = scalar[0]
            scalar_y = scalar[1]
            scalar_uy = scalar[2]
            
            
            
#            j2_x = j_2_all[0]
#            j2_y = j_2_all[1]
#            j2 =   j_2_all[2]
#            
            x = ht[0]
            y = ht[1]
            ht_array =ht[2] 
            
            # Points for triangulation
            points = [x,y]
            points = np.asarray(points)
            points = points.transpose()
            
            
            
            #------------------------------------------------------------------
            cmap_ht = 'RdBu'
            cmap_j2 = 'jet'
            cmap_uy = "RdBu"
            
            ax1 = plt.Subplot(fig,gs00[0,:])
            im = ax1.tripcolor(x,y,ht_array,shading='gouraud',cmap = cmap_ht,vmin = -15,vmax=15)
            
            # label for x and y axis
            ax1.set_xlabel('Distance from grounding line [m]')
            ax1.set_ylabel('Across flow direction [m]')
            start,end = ax1.get_xlim()
            
            # x-axis labelling for grounding-line distance
            ax1.xaxis.tick_bottom()
            ax1.xaxis.set_tick_params(pad=0)
            ax1.set_xticks([cs1,cs2,cs3])
            ax1.set_xticklabels([cs1-GL,cs2-GL,cs3-GL])
            
            # visibility of spines
            ax1.spines['top'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['left'].set_visible(False)
            #ax1.tick_params(direction='in',length=3,width=1)
         
            # Crosssections
            line_cs1 = ax1.vlines(x=cs1, ymin=ymin, ymax=ymax, color='r',linewidth=0.5)
            line_cs2 = ax1.vlines(x=cs2, ymin=ymin, ymax=ymax, color='r',linewidth=0.5)
            line_cs3 = ax1.vlines(x=cs3, ymin=ymin, ymax=ymax, color='r',linewidth=0.5)

            # Annotation for cross sections
            ax1.text(cs1,1800,"A",ha = 'center', va='top')
            ax1.text(cs1,-1800,"A'",ha = 'center', va='bottom')
            ax1.text(cs2,1800,"B",ha = 'center', va='top')
            ax1.text(cs2,-1800,"B'",ha = 'center', va='bottom')
            ax1.text(cs3,1800,"C",ha = 'center', va='top')
            ax1.text(cs3,-1800,"C'",ha = 'center', va='bottom')
            
            plt.setp(ax1.get_xticklabels(),fontweight = 'bold')
            plt.setp(ax1.get_yticklabels(),fontweight = 'bold')
            
            # Delaunay triangulation for grid
            delaunay = Delaunay(points)
            ax1.triplot(x,y,delaunay.simplices,alpha=0.05,linewidth=0.3)
            ax1.set_yticks([-2000,-1000,0,1000,2000])
            ax1.set_ylim(-2000,2000)
            ax1.set_xlim(1060000,1079000)
            plt.subplots_adjust(hspace=0.42)
            
            
            # place text box in upper left in axes coords      
            props = dict(boxstyle='round', facecolor='wheat')
            ax1.text(0.2, 0.92,'t = ' + str(t1) + 'a\n @ cw = ' + str(self.width*2) + 'm', transform=ax1.transAxes, 
                     verticalalignment='top', bbox=props,weight='bold')
            
            
            fig.add_subplot(ax1)
            
         
            #------------------------------------------------------------------
            ax2 = plt.Subplot(fig,gs00[1,0])
            cutted = mr.cut_and_slice(cs1,' syz')

            upper = c1[0]
            lower_2 = c1[1]
            
            original_2 = c1[2]
            ax2.plot(lower_2,'k-')
            
            #ax2.plot(upper,'b')
            ax2.plot(original_2,linestyle = '--', color = blue)
            ax2.pcolormesh(cutted[2],cutted[3],cutted[4],shading='gouraud',cmap = "bwr")
            print(np.max(cutted[4]))
#            legend_ax2 = ax2.legend(['Modelled thickness','Hydrostatic thickness'],loc="upper right", prop=dict(weight='bold'))
#            frame_ax2 = legend_ax2.get_frame()
#            frame_ax2.set_facecolor('0.7')
#            frame_ax2.set_edgecolor('0.7')
            
            ax2.set_ylabel('Shelf elevation [m]')
            ax2.set_xlim(-2000,2000)
            ax2.set_ylim([-350,-50])
            ax2.set_yticks([-300,-200,-100])
            bottom, top = ax2.get_ylim()
            ax2.text(ymin_2,top,"A",ha = 'right', va='bottom')
            ax2.text(ymax_2,top,"A'",ha = 'left', va='bottom')
            ax2.set_xticks([-1500,0,1500])
            fig.add_subplot(ax2)

            
            #------------------------------------------------------------------
            ax3 = plt.Subplot(fig,gs00[1,1])
            cutted = mr.cut_and_slice(cs2,' syz')
            
            
            upper = c2[0]
            lower_3 = c2[1]
            original_3 = c2[2]
            ax3.plot(lower_3,'k-')
            #ax3.plot (upper,'b')
            ax3.plot(original_3,linestyle = '--', color = blue)
           
            im_3 = ax3.pcolormesh(cutted[2],cutted[3],cutted[4],shading='gouraud',cmap = "bwr")
            print(np.max(cutted[4]))
                       
            ax3.set_xlim(-2000,2000)
            ax3.set_ylim([-350,-50])
            ax3.set_yticks([-250,-150,-50])
            bottom, top = ax3.get_ylim()
            ax3.set_xlabel('Across flow direction [m]')
            ax3.text(ymin_2,top,"B",ha = 'right', va='bottom')
            ax3.text(ymax_2,top,"B'",ha = 'left', va='bottom')
            ax3.axes.get_yaxis().set_visible(False)
            ax3.set_xticks([-1500,0,1500])
            fig.add_subplot(ax3)
            
            #------------------------------------------------------------------
            ax4 = plt.Subplot(fig,gs00[1,2])
            cutted = mr.cut_and_slice(cs3,' syz')
           
            
            upper = c3[0]
            lower_4 = c3[1]
            original_4 = c3[2]
            ax4.plot(lower_4,'k-')
            #ax4.plot (upper,'b')
            ax4.plot(original_4,linestyle = '--', color = blue)
            ax4.pcolormesh(cutted[2],cutted[3],cutted[4],shading='gouraud',cmap = "bwr")
            print(np.max(cutted[4]))
            
            ax4.set_xlim(-2000,2000)
            ax4.set_ylim([-350,-50])
            ax3.set_yticks([-250,-150,-50])
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
            ax4.axvline(0,linestyle='--',linewidth=0.5,color='grey')
            ax4_twin.axhline(0,linestyle='--',linewidth=0.5,color='grey')
                        
            ax3_twin = ax3.twinx()
            ax3_twin.plot(original_3-lower_3,color=red,linestyle=':')
            ax3_twin.set_yticks([-10,-5,0,5,10])
            ax3_twin.set_yticklabels([])
            ax3_twin.set_ylim(-15,15)
            ax3.axvline(0,linestyle='--',linewidth=0.5,color='grey')
            ax3_twin.axhline(0,linestyle='--',linewidth=0.5,color='grey')
            
            ax2_twin = ax2.twinx()
            ax2_twin.plot(original_2-lower_2,color=red,linestyle=':')
            ax2_twin.set_yticks([-10,-5,0,5,10])
            ax2_twin.set_yticklabels([])
            ax2_twin.set_ylim(-15,15)
            ax2.axvline(0,linestyle='--',linewidth=0.5,color='grey')
            ax2_twin.axhline(0,linestyle='--',linewidth=0.5,color='grey')
             
            plt.setp(ax4_twin.get_yticklabels(),color=red)
        
        
        #2D colorbars
        axin_2 = inset_axes(ax2,width='80%',height='10%',loc=9,bbox_to_anchor=(0,0.23,1,1),bbox_transform=ax2.transAxes)    
        cbar_2 = Colorbar(ax=axin_2, mappable = im_3, extend = 'both', \
            orientation = 'horizontal', ticklocation = 'top',ticks=[-0.005,0.005])   
        cbar_2.ax.tick_params(direction='in',length=1.5,width=0.5,pad=0.2)
        
        
        axin_3 = inset_axes(ax3,width='80%',height='10%',loc=9,bbox_to_anchor=(0,0.23,1,1),bbox_transform=ax3.transAxes)    
        cbar_3 = Colorbar(ax=axin_3, mappable = im_3, extend = 'both', \
            orientation = 'horizontal', ticklocation = 'top',ticks=[-0.005,0.005])           
        cbar_3.ax.tick_params(direction='in',length=1.5,width=0.5,pad=0.2)
        
        axin_4 = inset_axes(ax4,width='80%',height='10%',loc=9,bbox_to_anchor=(0,0.23,1,1),bbox_transform=ax4.transAxes)    
        cbar_4 = Colorbar(ax=axin_4, mappable = im_3, extend = 'both', \
            orientation = 'horizontal', ticklocation = 'top',ticks=[-0.005,0.005])          
        cbar_4.ax.tick_params(direction='in',length=1.5,width=0.5,pad=0.2)
        
        
        
        
        
        cbar_4.set_label('$\sigma_{yz}$ [MPa]',labelpad=0)
        cbar_3.set_label('$\sigma_{yz}$ [MPa]',labelpad=0)
        cbar_2.set_label('$\sigma_{yz}$ [MPa]',labelpad=0)
        
        
        #3D colorbar
        title_uy = "Velocity $u_{y}$ [$m\: a^{-1}$]"
        title_ht = 'Hydrostatic dev. [m]'
        title_j2 = 'Effective stress [MPa]'
        
        axin = inset_axes(ax1,width='100%',height='10%',loc=9,bbox_to_anchor=(0,0.15,1,1),bbox_transform=ax1.transAxes)    
        cbar = Colorbar(ax=axin, mappable = im, extend = 'both', \
        orientation = 'horizontal', ticklocation = 'top',ticks=MaxNLocator(3)) 
        cbar.set_label(title_ht,labelpad=0.1)
        cbar.ax.tick_params(direction='in',length=1,width=1,pad=0)
        
        
        
        
        
        
        
        path = str('plots/03_results/07_3D/')
        fname= str('3D_bridging_hd_' + str(self.amp) + str(self.width)+ '_' + str(t1) + '.png')
        fname_pdf = str('3d_bridging_hd_' + str(self.amp) + str(self.width)+ '_' + str(t1) + '.pdf')
        
        fig.savefig(path + fname, format = 'png', dpi=1000)
        fig.savefig(path + fname_pdf, format = 'pdf', dpi=1000)
                  
        plt.show()
        