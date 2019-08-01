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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

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
        
        # Custom params load from __plot_params
        plt.rcParams.update(params_bridging_2d) 
      
       
        
    
        # Coordinates for zoomed box
        x1,x2,y1,y2 = -2500,2500,-300,-200
        xx1,xx2,yy1,yy2 = -2500,2500,0,50
        
             
        self.original_width = int(round(np.sqrt(self.width)*2)) 
        nums = ["(a)","(b)","(c)","(d)"]
        
        
        fig = plt.figure()    
        
        gs_0 = gridspec.GridSpec(2,2, hspace = 0.08, wspace = 0.05, figure = fig)
        self.times = [t1, t2, t3, t4]
        
        
        self.deviation_list = []
        
        # 2,2 plots with 3 subplots
        for i ,t, num in zip(range(4), self.times, nums):
            gs00 = gridspec.GridSpecFromSubplotSpec(2,1, \
            height_ratios=[0.05,1], hspace = 0.02 ,subplot_spec=gs_0[i])
            
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
                cmap = str("bwr")
                title = str("Bridging $\sigma_{xy}$ [$MPa$]")
                scalarname = str('bridging')
            elif scalar is not None:
                self.scalar=scalar
                sxy = mr.get_scalar(str(self.scalar))
                sxy = abs(sxy[:,1])
                cmap = str("RdYlBu_r")
                title = str("Velocity $u_{y}$ [$m\: a^{-1}$]")
                scalarname = str('velo')
                
            points = ht[4]
            self.x = points[:,0]
            self.y = points[:,1]
            
            x = points[:,0]
            y = points[:,1]
            
         
            new_points = x,y
    
            self.xx,self.yy = np.meshgrid(x,y)
       
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
           
            im = ax1.pcolormesh(self.xx,self.yy,grid,vmin=0,vmax=0.8,cmap = cmap)
            
            # Set labels
            ax1.set_xlabel('Distance [m]', visible = False)
            ax1.set_ylabel('Height [m]',visible=False)
            #ax1.get_yaxis().set_ticklabels([])
            #ax1.get_xaxis().set_ticklabels([])
            ax1.set_xticks([-1000,0,1000])
            ax1.set_yticks([-250,-100,0,50])
            
            plt.setp(ax1.get_xticklabels(),fontweight = 'bold')
            plt.setp(ax1.get_yticklabels(),fontweight = 'bold')
            
            upper=ht[2]
            lower=ht[3]
            points=ht[4]
            x_line = ht[0]
            
            # Plot lines
            ax1.plot(x_line,lower,'b-',linewidth=1)
            ax1.fill_between(x_line, lower, lower.min(), color='w')
            ax1.plot(x_line,upper,'b-',linewidth=1)
            ax1.fill_between(x_line, upper, upper.max(), color='w')
            ax1.tick_params(direction='in',length=3,width=1)
            ax1.set_xlim(-2500,2500)
            ax1.set_ylim(-290,100)
            
            # Hydrostatic thickness
            ax1.plot(self.new_x,self.calc_thickness_bs,'r--',linewidth=1)
            
            # Legend
            legend_ax2 = ax1.legend(['Modelled thickness for t=' + str(t*5)+'a','_nolegend_','Hydrostatic thickness'],loc="lower left",bbox_to_anchor=(0.28,0.82))
            frame_ax2 = legend_ax2.get_frame()
            frame_ax2.set_facecolor('0.7')
            frame_ax2.set_edgecolor('0.7')
            
            # Colorbar
            cbar = plt.subplot(gs00[0,:])
            cbar = Colorbar(ax=cbar, mappable = im, extend = 'both', \
            orientation = 'horizontal', ticklocation = 'top')
            cbar.ax.xaxis.set_tick_params(pad=1,size=1.5)
           
            cbar.ax.xaxis.set_ticks([sxy.max()])
            cbar.ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
            
            # place text box in upper left in axes coords
            
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
            
            ax1.text(0.05, 0.95, num, transform=ax1.transAxes, fontsize=7,
                verticalalignment='top', bbox=props)
            
            
            
            
            if i==0:
                ax1.get_xaxis().set_ticklabels([]) 
            if i==1:
                ax1.get_yaxis().set_ticklabels([]) 
                ax1.get_xaxis().set_ticklabels([]) 
                
               
            if i==1:    
                ax1.get_xaxis().set_ticklabels([]) 
            
            if i==3:    
                ax1.get_yaxis().set_ticklabels([]) 
                
            
            
            
            
            plt.suptitle(title + ' @ cw = ' + str(self.original_width*2) + 'm', weight='bold', fontsize = 7.5, y = 0.95)
            fig.text(0.5,0.06, 'Channel-width [m]',ha = 'center',fontsize=8)
            fig.text(0.05,0.6, 'Shelf-height [m]',ha = 'center',rotation = 'vertical',fontsize=8)
            fig.add_subplot(ax1)
        
                 
            
#            #------------------------------------------------------------------         
#            # Plot calculated (red) and modelled (blue) hydrostatic thickness 
#            #------------------------------------------------------------------
#            ax2 = plt.Subplot(fig,gs00[2,:])
#            
#            ax2.plot(self.new_x,self.lower,'b-', label = 'Modelled thickness',linewidth=1)
#            ax2.plot(self.new_x,self.upper,'b',linewidth=1)
#            ax2.plot(self.new_x,self.calc_thickness_bs,'r--',label = 'Hydrostatic thickness',linewidth=1)
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
#            ax2.tick_params(direction='in',length=3,width=1)
            
           
            
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
       
        
       
        
        path = str('plots/')
        fname_png = str('_dev_2d__' + str(self.original_width) + '.png')
       # fname_pdf = str('briding_dev_2d__' + str(self.original_width) + '.pdf')
    
        plt.savefig(path + scalarname + fname_png, format = 'png',dpi=1000, bbox_inches = 'tight')
        #plt.savefig(path + fname_pdf, format = 'pdf',dpi=1000, bbox_inches = 'tight')
        
        plt.show() 
        
  