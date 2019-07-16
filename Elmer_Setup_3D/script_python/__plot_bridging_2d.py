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

from __main__ import ModelRun
from __plot_params import params


from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
import matplotlib.font_manager

import numpy as np
import matplotlib.pyplot as plt
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
                 scalar = None,
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
        # Setup input
        #======================================================================
        
        self.width = width
            
     
        
        #======================================================================
        # Define variables for input
        #======================================================================
       
              
        # Choose default timesteps if timesteps are not given
        if (t1,t2,t3,t4) is None:
                t1 = 50
                t2 = 100
                t3 = 150
                t4 = 200
        else:
                self.t1 = t1
                self.t2 = t2
                self.t3 = t2
                self.t4 = t3        
                        
            
            
        if prop is None:
            prop = 0        
        else:
            self.prop = prop  
        
        
        # Boundaries ax2 plot
        ymin = -350
        ymax = 60
        
      
        
        #======================================================================
        # Define figure and plot properties
        #======================================================================
        # Custom model load from __plot_params
        plt.rcParams.update(params) 
      
       
        
    
        # Coordinates for zoomed box
        x1,x2,y1,y2 = -2500,2500,-300,-200
        xx1,xx2,yy1,yy2 = -2500,2500,0,50
        
             
        self.original_width = int(round(np.sqrt(self.width)*2)*2) 
        
         
        fig = plt.figure(figsize = (40,50))    
        
        gs_0 = gridspec.GridSpec(2,2, hspace = 0.2, wspace = 0.2, figure = fig)
        self.times = [t1, t2, t3, t4]
        
        
        self.deviation_list = []
        
        # 2,2 plots with 3 subplots
        for i ,t in zip(range(4), self.times):
            gs00 = gridspec.GridSpecFromSubplotSpec(3,3, \
            height_ratios=[0.05,1,0], hspace = 0.10 ,subplot_spec=gs_0[i])
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
            
            if scalar is None:
                sxy = mr.sxy            
            else:
                self.scalar=scalar
                sxy = mr.get_scalar(str(self.scalar))
                sxy = sxy
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
            #ax1(frameon=0)
            im = ax1.pcolormesh(self.xx,self.yy,grid,vmin=-0.03,vmax=0.03,cmap = 'RdBu')
            
            # Colorbar definition (extra axis)
            cbar = plt.subplot(gs00[0,:])
            cbar = Colorbar(ax=cbar, mappable = im, extend = 'both', \
            orientation = 'horizontal', ticklocation = 'top') 
            cbar.set_label('Bridging stress $\sigma_{xy}$ \n' \
            + 't = ' + str(t*5)+ ' a' + ', channel width = ' + str(self.original_width) + 'm \n' + \
            '$\sigma_{xy} max$'+ ' = ' + str(max_sxy.round(4)) + ' MPa',labelpad=8)
            
            ax1.set_xlabel('Distance [m]')
            ax1.set_ylabel('Height [m]')
            ax1.minorticks_on()
            upper=ht[2]
            lower=ht[3]
            points=ht[4]
            x_line = ht[0]
            ax1.plot(x_line,lower,'b-',linewidth=5)
            ax1.plot(x_line,upper,'b-',linewidth=5)
            ax1.tick_params(direction='in',length=6,width=2)
            ax1.set_xlim(-2500,2500)
            
            fig.add_subplot(ax1)
            
            
            
#            #------------------------------------------------------------------         
#            # Plot calculated (red) and modelled (blue) hydrostatic thickness 
#            #------------------------------------------------------------------
#            ax2 = plt.Subplot(fig,gs00[2,:])
#            
#            ax2.plot(self.new_x,self.lower,'b-', label = 'Modelled thickness',linewidth=5)
#            ax2.plot(self.new_x,self.upper,'b',linewidth=5)
#            ax2.plot(self.new_x,self.calc_thickness_bs,'r--',label = 'Hydrostatic thickness',linewidth=5)
#           
#            # Legend
#            legend = ax2.legend(loc = 'lower right', bbox_to_anchor=(0.3,0.76))
#            frame = legend.get_frame()
#            frame.set_facecolor('0.9')
#            frame.set_edgecolor('0.9')
#            
#            ax2.grid()
#            ax2.set_xlabel('Distance [m]')
#            ax2.set_ylabel('Height [m]')
#            ax2.minorticks_on()
#            ax2.set_ylim([ymin,ymax])
#            bottom, top = ax2.get_ylim()
#            ax2.tick_params(direction='in',length=6,width=2)
#            
#           
#            
#            fig.add_subplot(ax2)
#            
#            
#            # Zoomed in rectangle for better visualisation of channel
#            axins = zoomed_inset_axes(ax2,1.5,loc='upper right',borderpad=5)
#            axins.plot(self.new_x,self.calc_thickness_bs,'r--',linewidth=5)
#            axins.plot(self.new_x,lower,'b-',linewidth=5)
#            axins.set_xlim(x1,x2)
#            axins.set_ylim(y1,y2)
#            axins.minorticks_on()
#            axins.grid(which='minor',linestyle = ':', linewidth = '0.5', color = 'black')
#            
#            # Visibility of ticks
#            #plt.yticks(visible=False)
#            #plt.xticks(visible=False)
#            
#            # Connection to bounding box
#           #mark_inset(ax2,axins,loc1=3,loc2=4,fc="none", ec="0.5")
#            
#            
#            # Zoomed in rectangle for better visualisation of surface expression
#            axins = zoomed_inset_axes(ax2,1.5,loc='upper left',borderpad=5)            
#            axins.plot(self.new_x,upper,'b-')
#            axins.set_xlim(xx1,xx2)
#            axins.set_ylim(yy1,yy2)
#            axins.minorticks_on()
#            axins.grid(which='minor',linestyle = ':', linewidth = '0.5', color = 'black')
#            
#            # Visibility of ticks
#            #plt.yticks(visible=False)
#            #plt.xticks(visible=False)
#            
#            # Connection to bounding box
#            #mark_inset(ax2,axins,loc1=3,loc2=4,fc="none", ec="0.5")
           
    
                  
    
    
#        path = str('plots/')
#        fname= str('briding_dev_2d__' + str(original_width) + '.eps')
#        
#        fig.savefig(path + fname, format = 'eps',dpi=1000)
        self.fig = fig
        fig.show()
            
    def getfig(self):
        path = str('plots/')
        fname= str('briding_dev_2d__' + str(self.original_width) + '.eps')
    
        self.fig.savefig(path + fname, format = 'eps',dpi=1000)