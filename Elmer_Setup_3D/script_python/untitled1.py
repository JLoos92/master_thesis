#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 16:01:02 2019

@author: jloos
"""

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

from main import ModelRun
from __plot_params import params_bridging_2d


from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
import matplotlib.font_manager

import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy.interpolate import griddata,interpolate






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
        plt.rcParams.update(params_bridging_2d) 
      
       
        
    
        # Coordinates for zoomed box
        x1,x2,y1,y2 = -2500,2500,-300,-200
        xx1,xx2,yy1,yy2 = -2500,2500,0,50
        
             
        self.original_width = int(round(np.sqrt(self.width)*2)) 
        
         
        fig = plt.figure(frameon=False)    
        
        gs_0 = gridspec.GridSpec(2,2, hspace = 0.05, wspace = 0.05, figure = fig)
        self.times = [t1, t2, t3, t4]
        
        subs = []
        self.deviation_list = []
        
        # 2,2 plots with 3 subplots
        for i ,t in zip(range(4), self.times):
            gs00 = gridspec.GridSpecFromSubplotSpec(3,3, \
            height_ratios=[0.05,1,0.2], hspace = 0.2 ,subplot_spec=gs_0[i])
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
    
            self.xx,self.yy = np.meshgrid(200,200)
       
            grid = griddata(new_points,sxy,(self.xx,self.yy),method='cubic')
            #print(grid)
           
            
            # Points for triangulation
            points = [self.x,self.y]
            points = np.asarray(points)
            points = points.transpose()
            max_sxy = np.max(sxy)
            min_sxy = np.min(sxy)
            
            
      
            #------------------------------------------------------------------
            ax1 = plt.Subplot(fig,gs00[1,:])
            #ax1(frameon=0)
            im = ax1.pcolormesh(self.xx,self.yy,grid,vmin=-0.03,vmax=0.03,cmap = 'RdBu')
            
#            # Colorbar definition (extra axis)
#            cbar = plt.subplot(gs00[0,:])
#            cbar = Colorbar(ax=cbar, mappable = im, extend = 'both', \
#            orientation = 'horizontal', ticklocation = 'top') 
#            cbar.set_label('Bridging stress $\sigma_{xy}$ \n' \
#            + 't = ' + str(t*5)+ ' a' + ', channel width = ' + str(self.original_width) + 'm \n' + \
#            '$\sigma_{xy} max$'+ ' = ' + str(max_sxy.round(4)) + ' MPa',labelpad=8)
            
            ax1.set_xlabel('Distance [m]', visible = False)
            ax1.set_ylabel('Height [m]',visible=False)
            ax1.get_yaxis().set_ticklabels([])
            ax1.get_xaxis().set_ticklabels([])
            ax1.set_xticks([-2000,-self.original_width,0,self.original_width,2000])
            ax1.set_yticks([-250,-100,0,50])
            ax1.minorticks_on()
            upper=ht[2]
            lower=ht[3]
            points=ht[4]
            x_line = ht[0]
            
            # Plot lines
            ax1.plot(x_line,lower,'b-',linewidth=2)
            ax1.fill_between(x_line, lower, lower.min(), color='w')
            ax1.plot(x_line,upper,'b-',linewidth=2)
            ax1.tick_params(direction='in',length=3,width=1)
            ax1.set_xlim(-2500,2500)
            labels = [-2000,-1000,0,1000,2000]
            
           
            
            
            fig.add_subplot(ax1)
            
            

        
        
        path = str('plots/bridging_2d')
        fname_png = str('briding_dev_2d__' + str(self.original_width) + '.png')
       
    
        plt.savefig(path + fname_png, format = 'png',dpi=1000,bbox_inches = 'tight')
        
        plt.show() 
        
