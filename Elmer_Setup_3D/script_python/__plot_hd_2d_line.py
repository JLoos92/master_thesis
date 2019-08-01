#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:29:50 2019

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

class Plot_line_hd_2d():
    
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
        
        #======================================================================
        # Define figure and plot properties
        #======================================================================
        
        # Custom params load from __plot_params
        plt.rcParams.update(params_bridging_2d) 
          
        self.original_width = int(round(np.sqrt(self.width)*2)) 
        nums = ["(a)","(b)","(c)","(d)"]
        
        
        fig = plt.figure()    
        
        gs_0 = gridspec.GridSpec(2,2, hspace = 0.01, wspace = 0.05, figure = fig)
        self.times = [t1, t2, t3, t4]
        
        
        self.deviation_list = []
        
        # 2,2 plots with 3 subplots
        for i ,t, num in zip(range(4), self.times, nums):
            gs00 = gridspec.GridSpecFromSubplotSpec(2,1, \
            height_ratios=[0.05,1], hspace = 0.01 ,subplot_spec=gs_0[i])
            
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
              
            points = ht[4]
            self.x = points[:,0]
            self.y = points[:,1]
            
            x = points[:,0]
            y = points[:,1]
            
         
            new_points = x,y
           
            # Points for triangulation
            points = [self.x,self.y]
            points = np.asarray(points)
            points = points.transpose()

            ax1 = plt.Subplot(fig,gs00[1,:])

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
            
            
            
            # compute deviation regular
            lower = ht[3]
            calc_thickness_bs = ht[1]


            
            hydrostatic_deviation = calc_thickness_bs - lower
            
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
                
            
            
            
            title = str('Hydrostatic deviation')
            plt.suptitle(title + ' @ cw = ' + str(self.original_width*2) + 'm', weight='bold', fontsize = 7.5, y = 0.90)
            fig.text(0.5,0.06, 'Channel-width [m]',ha = 'center',fontsize=8)
            fig.text(0.05,0.6, 'Shelf-height [m]',ha = 'center',rotation = 'vertical',fontsize=8)
            fig.add_subplot(ax1)
            
            
            ax2 = ax1.twinx()
            ax2.plot(x_line, hydrostatic_deviation, 'g:',linewidth=1.5)
            
            if i==0:
                ax2.get_yaxis().set_ticklabels([])
            if i==2:    
                ax2.get_yaxis().set_ticklabels([])
                  
            ax2.get_yaxis().set_ticks([-10,-5,0,5,10])
            ax2.set_ylim(-30,30)
            plt.setp(ax2.get_yticklabels(),fontweight = 'bold',color='g')
            #ax2.set_ylabel('$\sigma_{xy}$ [MPa]',labelpad = 15, color = 'b')
            fig.text(1,0.6, 'Hydrostatic deviation [m]',ha = 'center',rotation = 'vertical',fontsize=8,color='g')
       
        
        path = str('plots/Final_plots/hd_line_2d/')
        fname_png = str('line_dev_2d__' + str(self.original_width) + '.png')
        fname_pdf = str('line_dev_2d__' + str(self.original_width) + '.pdf')
    
        plt.savefig(path + fname_png, format = 'png',dpi=1000, bbox_inches = 'tight')
        plt.savefig(path + fname_pdf, format = 'pdf',dpi=1000, bbox_inches = 'tight')
        
        plt.show() 
