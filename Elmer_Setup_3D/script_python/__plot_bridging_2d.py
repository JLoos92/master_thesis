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
        self.width = width
        
        # Boundaries for model domain
        ymin = -350
        ymax = 60

        
        original_width = round(np.sqrt(self.width)*2)*2
        
        
        #======================================================================
        # Define figure and plot properties
        #======================================================================
        
        fig = plt.figure(figsize = (40,40))    
        
        gs_0 = gridspec.GridSpec(2,2, hspace = 0.2, wspace = 0.1, figure = fig)
        self.times = [t1, t2, t3, t4]
        
        # 2,2 plots with 3 subplots
        for i ,t in zip(range(4), self.times):
            gs00 = gridspec.GridSpecFromSubplotSpec(3,3, \
            height_ratios=[0.05,1,0.8], hspace = 0.25 ,subplot_spec=gs_0[i])
        
        # Run ModelRun-Class with default or self
            if None is (self.width):
                mr = ModelRun(100000,0,0,0,t)
            mr = ModelRun(self.width,0,0,0,t,"vtu")
            
           
            ht = mr.compute_hydrostatic_thickness()
            calc_thickness_bs = ht[1]
            
            points = ht[4]
            x = points[:,0]
            y = points[:,1]
           
            #ht_array =ht[2] 
            sxy = mr.sxy
            # Points for triangulation
            points = [x,y]
            points = np.asarray(points)
            points = points.transpose()
            max_sxy = np.max(sxy)
            min_sxy = np.min(sxy)
            masked_sxy = ma.masked_where(sxy>-10,sxy)
            print(masked_sxy)
            
            
            #------------------------------------------------------------------
            ax1 = plt.Subplot(fig,gs00[1,:])
            im = ax1.tripcolor(x,y,masked_sxy,shading="gouraud",vmin=-0.03,vmax=0.03,cmap = 'RdBu')
            
            # Colorbar definition (extra axis)
            cbar = plt.subplot(gs00[0,:])
            cbar = Colorbar(ax=cbar, mappable = im, extend = 'both', \
            orientation = 'horizontal', ticklocation = 'top') 
            cbar.set_label('Bridging stress $\sigma_{xy}$ \n' \
            + 't = ' + str(t*10) + ', channel width = ' + str(original_width) + 'm \n' + \
            '$\sigma_{xy} max$'+ ' = ' + str(max_sxy.round(4)),labelpad=10, fontdict = font_title)
            
            ax1.set_xlabel('Width [m]',fontdict = font_label)
            ax1.set_ylabel('Height [m]',fontdict = font_label)
            ax1.minorticks_on()
            upper=ht[2]
            lower=ht[3]
            points=ht[4]
            x_line = ht[0]
            ax1.plot(x_line,lower,'b-')
            ax1.plot(x_line,upper,'r-')


         
            # Delaunay triangulation for grid
            #delaunay = Delaunay(points)
            #ax1.triplot(x,y,delaunay.simplices,alpha=0.05)
            
            fig.add_subplot(ax1)
            
#            # Zoomed in rectangle for better visualisation of channel
#            axins = zoomed_inset_axes(ax1,1.5,loc='upper right',borderpad=5)
#            axins.tripcolor(x,y,ht_array,shading='gouraud',vmin=-15,vmax=15,cmap = 'RdBu')
#            axins.set_xlim(x1,x2)
#            axins.set_ylim(y1,y2)
#            plt.yticks(visible=False)
#            plt.xticks(visible=False)
#            mark_inset(ax1,axins,loc1=2,loc2=4,fc="none", ec="0.5")
                
           
            
            #------------------------------------------------------------------
            ax2 = plt.Subplot(fig,gs00[2,:])
            upper=ht[2]
            lower=ht[3]
            #points=ht[4]
            x = ht[0]
            ax2.plot(x,lower,'b-', label = 'Modelled thickness')
            ax2.plot(x,upper,'b')
            ax2.plot(x,calc_thickness_bs,'r--',label = 'Hydrostatic thickness')
           # ax2.plot(original,'r--', label = 'Hydrostatic thickness')
            ax2.legend(loc = 'lower right', bbox_to_anchor=(0.55,0.77))
            ax2.grid()
            ax2.set_xlabel('Width [m]', fontdict = font_label)
            ax2.set_ylabel('Height [m]', fontdict = font_label)
            ax2.minorticks_on()
            ax2.set_ylim([ymin,ymax])
            bottom, top = ax2.get_ylim()
            ax2.text(ymin,top,"A",fontdict=font_annotation,ha = 'right', va='bottom')
            ax2.text(ymax,top,"A'",fontdict=font_annotation,ha = 'left', va='bottom')
            fig.add_subplot(ax2)
           

#            #------------------------------------------------------------------
#            ax3 = plt.Subplot(fig,gs00[2,1])
#            lower = c2[0]
#            upper = c2[1]
#            original = c2[2]
#            ax3.plot(lower,'b-')
#            ax3.plot (upper,'b')
#            ax3.plot(original,'r--')
#            
#            ax3.grid()
#            ax3.set_ylim([-350,55])
#            bottom, top = ax3.get_ylim()
#            ax3.text(ymin,top,"B",fontdict=font_annotation,ha = 'right', va='bottom')
#            ax3.text(ymax,top,"B'",fontdict=font_annotation,ha = 'left', va='bottom')
#            fig.add_subplot(ax3)
#            
#            #------------------------------------------------------------------
#            ax4 = plt.Subplot(fig,gs00[2,2])
#            lower = c3[0]
#            upper = c3[1]
#            original = c3[2]
#            ax4.plot(lower,'b-')
#            ax4.plot (upper,'b')
#            ax4.plot(original,'r--')
#            
#            ax4.grid()
#            ax4.set_ylim([-350,55])
#            bottom, top = ax4.get_ylim()
#            ax4.text(ymin, top,"C",fontdict=font_annotation,ha = 'right', va='bottom')
#            ax4.text(ymax, top,"C'",fontdict=font_annotation,ha = 'left', va='bottom')
#            fig.add_subplot(ax4)
            
        
                  
        fig.show()
        