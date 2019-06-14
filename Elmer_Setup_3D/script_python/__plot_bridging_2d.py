#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:05:34 2019

@author: jloos
"""

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
import numpy.ma as ma
from scipy.interpolate import griddata






class Plot_bridging_2d():
    
    """
    Class Plot_hydrostatic_deviation:
    Here description
    """
    

    
    
    
    def __init__(self,
                 t1 = None,
                 t2 = None,
                 t3 = None,
                 t4 = None,
                 width = None,
                 prop = None,
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
        
        if prop is None:
            prop = 0        
                
        self.t1 = t1
        self.t2 = t2
        self.t3 = t2
        self.t4 = t3
        self.width = width
        self.prop = prop
        
        # Boundaries ax2 plot
        ymin = -350
        ymax = 60
        
      
        
        #======================================================================
        # Define figure and plot properties
        #======================================================================
       
      
       
        

        # Coordinates for zoomed box
        x1,x2,y1,y2 = -2500,2500,-300,-200
        xx1,xx2,yy1,yy2 = -2500,2500,0,50
        
             
        original_width = round(np.sqrt(self.width)*2)*2
        
         
        fig = plt.figure(figsize = (40,40))    
        
        gs_0 = gridspec.GridSpec(2,2, hspace = 0.2, wspace = 0.1, figure = fig)
        self.times = [t1, t2, t3, t4]
        
        
        self.deviation_list = []
        
        # 2,2 plots with 3 subplots
        for i ,t in zip(range(4), self.times):
            gs00 = gridspec.GridSpecFromSubplotSpec(3,3, \
            height_ratios=[0.05,1,0.8], hspace = 0.25 ,subplot_spec=gs_0[i])
             # Run ModelRun-Class with default or self
            if None is (self.width):
                mr = ModelRun(100000,0,0,0,t)
            mr = ModelRun(150,self.width,0,self.prop,t,"2")
            
           
            ht = mr.compute_hydrostatic_thickness()
            
            self.upper=ht[2]
            self.lower=ht[3]
            self.new_x = ht[0]
        
            # Calculated hydrostatic thickness
            self.calc_thickness_bs = ht[1]
            
            sxy = mr.sxy
            points = ht[4]
            self.x = points[:,0]
            self.y = points[:,1]
            
            x = points[:,0]
            y = points[:,1]

            new_points = x,y

            self.xx,self.yy = np.meshgrid(x,y)

            
            grid = griddata(new_points,sxy,(self.xx,self.yy),method='nearest')
            
           
            
            # Points for triangulation
            points = [self.x,self.y]
            points = np.asarray(points)
            points = points.transpose()
            max_sxy = np.max(sxy)
            min_sxy = np.min(sxy)
            masked_sxy = ma.masked_where(sxy>-10,sxy)
            print(masked_sxy)
      
            #------------------------------------------------------------------
            ax1 = plt.Subplot(fig,gs00[1,:])
            ax1(frameon=0)
            im = ax1.pcolormesh(self.xx,self.yy,grid,vmin=-0.03,vmax=0.03,cmap = 'RdBu')
            
            # Colorbar definition (extra axis)
            cbar = plt.subplot(gs00[0,:])
            cbar = Colorbar(ax=cbar, mappable = im, extend = 'both', \
            orientation = 'horizontal', ticklocation = 'top') 
            cbar.set_label('Bridging stress $\sigma_{xy}$ \n' \
            + 't = ' + str(t*10) + ', channel width = ' + str(original_width) + 'm \n' + \
            '$\sigma_{xy} max$'+ ' = ' + str(max_sxy.round(4)),labelpad=10, fontdict = font_title)
            
            ax1.set_xlabel('Distance [m]',fontdict = font_label)
            ax1.set_ylabel('Height [m]',fontdict = font_label)
            ax1.minorticks_on()
            upper=ht[2]
            lower=ht[3]
            points=ht[4]
            x_line = ht[0]
            ax1.plot(x_line,lower,'b-')
            ax1.plot(x_line,upper,'b-')
            
          
            fig.add_subplot(ax1)
            
            
            
            #------------------------------------------------------------------         
            # Plot calculated (red) and modelled (blue) hydrostatic thickness 
            #------------------------------------------------------------------
            ax2 = plt.Subplot(fig,gs00[2,:])
            
            ax2.plot(self.new_x,self.lower,'b-', label = 'Modelled thickness')
            ax2.plot(self.new_x,self.upper,'b')
            ax2.plot(self.new_x,self.calc_thickness_bs,'r--',label = 'Hydrostatic thickness')
           
            # Legend
            legend = ax2.legend(loc = 'lower right', bbox_to_anchor=(0.55,0.77))
            frame = legend.get_frame()
            frame.set_facecolor('0.9')
            frame.set_edgecolor('0.9')
            
            ax2.grid()
            ax2.set_xlabel('Distance [m]', fontdict = font_label)
            ax2.set_ylabel('Height [m]', fontdict = font_label)
            ax2.minorticks_on()
            ax2.set_ylim([ymin,ymax])
            bottom, top = ax2.get_ylim()
            
            
           
            
            fig.add_subplot(ax2)
            
            
            # Zoomed in rectangle for better visualisation of channel
            axins = zoomed_inset_axes(ax2,1.5,loc='upper right',borderpad=5)
            axins.plot(self.new_x,self.calc_thickness_bs,'r--')
            axins.plot(self.new_x,lower,'b-')
            axins.set_xlim(x1,x2)
            axins.set_ylim(y1,y2)
            axins.minorticks_on()
            axins.grid(which='minor',linestyle = ':', linewidth = '0.5', color = 'black')
            
            # Visibility of ticks
            #plt.yticks(visible=False)
            #plt.xticks(visible=False)
            
            # Connection to bounding box
           #mark_inset(ax2,axins,loc1=3,loc2=4,fc="none", ec="0.5")
            
            
            # Zoomed in rectangle for better visualisation of surface expression
            axins = zoomed_inset_axes(ax2,1.5,loc='upper left',borderpad=5)            
            axins.plot(self.new_x,upper,'b-')
            axins.set_xlim(xx1,xx2)
            axins.set_ylim(yy1,yy2)
            axins.minorticks_on()
            axins.grid(which='minor',linestyle = ':', linewidth = '0.5', color = 'black')
            
            # Visibility of ticks
            #plt.yticks(visible=False)
            #plt.xticks(visible=False)
            
            # Connection to bounding box
            #mark_inset(ax2,axins,loc1=3,loc2=4,fc="none", ec="0.5")
           

                  
        fig.show()
        
        
        
        
#
            
            
            