#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:29:50 2019

@author: jloos
"""
from main import ModelRun
from __plot_params import params_horizontal


from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
import matplotlib.font_manager
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy.interpolate import griddata,interpolate

   

def plot_hd_line_2d(t1 = None,
                    t2 = None,
                    width = None,
                    prop = None):
                
    '''
    Function: plot_hd_line_2d
    ---------------
    Introduction figure for the deviation of hydrostatic equilibrium and bottom 
    + top surface. Two subplots with two timesteps for a given width. 
    Annotation and vertical lines for seperation are included.
    
    Parameters
    ----------
    t1 : int
        first timestep for plot
    t2 : int
        second timestep for plot
    width : int (w_mod)
        timestep for plot
    prop : int (0 for regular domain)
        timestep for plot    
    '''
       
    
    # Choose default timesteps if timesteps are not given    
    if (t1,t2) is None:
        t1 = 50
        t2 = 100
         
                    
    if prop is None:
        prop = 0        
      
    # Custom color definition    
    orange = '#D55E00'
    blue = '#396AB1'    
    red = '#CC2529'
    wheat = '#948B3D' 
    
    
    # Calculate original width
    original_halfwidth = int(round(np.sqrt(width)*2))
    original_width = original_halfwidth
    
    # Numeration and timesteps
    nums = ["(a)","(b)"]
    times = [t1, t2]
    

    # Set figure layout
    nrow = 2; ncol = 1;
    fig, ax1 = plt.subplots(nrows=nrow, ncols=ncol, sharex = True, squeeze = False, constrained_layout=True)

    # For common axis add subplot
    ax_common = fig.add_subplot(111,frameon=False)   
    
    
    plt.tick_params(labelcolor=False,top=False,bottom=False,left=False,right=False)

    # Custom params load from __plot_params
    plt.rcParams.update(params_horizontal) 
    
    
    for ax,num,t in zip(ax1.reshape(-1),nums,times): 
        
        
            mr = ModelRun(150,width,0,prop,t,"2")
            ht = mr.compute_hydrostatic_thickness()
            
            upper=ht[2]
            lower=ht[3]
            new_x = ht[0]
        
            # Calculated hydrostatic thickness
            calc_thickness_bs = ht[1]
              
            points = ht[4]
            x = points[:,0]
            y = points[:,1]

            # Points for triangulation
            points = [x,y]
            points = np.asarray(points)
            points = points.transpose()

            
            # Set labels
            ax.set_xticks([-original_width,0,original_width])
            ax.set_yticks([-250,-200,-100,-50,0,50])
            plt.setp(ax.get_xticklabels(),fontweight = 'bold')
            plt.setp(ax.get_yticklabels(),fontweight = 'bold')
            
            upper=ht[2]
            lower=ht[3]
            points=ht[4]
            x_line = ht[0]
            
            # Plot lines
            ax.plot(x_line,lower,'k-',linewidth=1.5)
            ax.plot(x_line,upper,'k-',linewidth=1.5,label='_nolegend_')
            ax.fill_between(x_line, lower, lower.min(), color='w')

            ax.fill_between(x_line, upper, upper.max(), color='w')
            ax.tick_params(direction='in',length=4,width=2)
            ax.set_xlim(-1500,1500)
            ax.set_ylim(-290,50)
            
            

            
            # Hydrostatic thickness
            ax.plot(new_x,calc_thickness_bs, linestyle = "--", color = blue ,linewidth=1.5)
            
            # Legend
            legend_ax2 = ax.legend(['Modelled thickness','Hydrostatic thickness'],loc="center right", prop=dict(weight='bold'))
            frame_ax2 = legend_ax2.get_frame()
            frame_ax2.set_facecolor('0.7')
            frame_ax2.set_edgecolor('0.7')
                        
            
            # compute deviation regular
            lower = ht[3]
            calc_thickness_bs = ht[1]
            
            hydrostatic_deviation = calc_thickness_bs - lower
            
            # place text box in upper left in axes coords      
            props = dict(boxstyle='round', facecolor='wheat')
          
            # Second axis
            axs = ax.twinx()
            axs.tick_params(direction='in',length=4,width=2)
              
            axs.set_yticks([-30,-10,0,10,30])
            axs.set_ylim(-40,40)
            
            plt.setp(axs.get_yticklabels(),fontweight = 'bold',color= red)
            
            axs.plot(x_line, hydrostatic_deviation, linestyle = ':', color = red)
            # Legend
            legend_ax2 = ax.legend(['Modelled thickness','Hydrostatic thickness'],bbox_to_anchor=[0.8,0.3],loc="center", prop=dict(weight='bold'))
            frame_ax2 = legend_ax2.get_frame()
            frame_ax2.set_facecolor('0.7')
            frame_ax2.set_edgecolor('0.7')


            title = str('Hydrostatic deviation')
            title_upper = str('Upper model-domain (zs)')
            
            ax.text(0.04, 0.92, num + '  t = ' + str(t*5) + 'a\n @ cw = ' + str(original_width*2) + 'm', transform=ax.transAxes, 
                verticalalignment='top', bbox=props, weight='bold')
            
            # Spines in subplot
            axs.spines['top'].set_visible(False)
            #axs.spines['bottom'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            
            fig.add_subplot(ax)
            
            

    # Set legend downer off
    ax1[1][0].get_legend().remove()    
    
    #Annotation
    ax1[1][0].annotate('Channel peak deviation',xy=(0,-242),xytext=(0,-50),va='bottom',ha='center',
                        arrowprops=dict(facecolor='black',shrink=0.01,width=0.5,headwidth=4,headlength=3.2))
    
    # Vertical lines
    ax1[1][0].vlines(-original_width,-290,-50,linewidth=0.5)
    ax1[1][0].vlines(original_width,-290,-50,linewidth=0.5)
    ax1[1][0].vlines(-original_width,-20,0,linewidth=0.5)
    ax1[1][0].vlines(original_width,-20,0,linewidth=0.5)
    
    ax1[0][0].vlines(-original_width,-290,-55,linewidth=0.5)
    ax1[0][0].vlines(original_width,-290,-55,linewidth=0.5)
    ax1[0][0].vlines(-original_width,-30,0,linewidth=0.5)
    ax1[0][0].vlines(original_width,-30,0,linewidth=0.5)
    
    
    
    # Annotations, labels etc.
    ax1[0][0].text(0,-120,'ICD',ha='center')
    ax1[1][0].text(-1000,-200,'OCD',ha='center')
    ax1[1][0].text(1000,-200,'OCD',ha='center')
    
    ax1[0][0].annotate('',xy=(-420,-85),xytext=(420,-85),
                        arrowprops=dict(facecolor='black',arrowstyle='<->'))
    ax1[0][0].text(0,-55,' Channel abutment deviation',ha='center')
    ax1[0][0].annotate('Outer peak deviation',xy=(-780,-168),xytext=(-780,-230),va='bottom',ha='center',
                        arrowprops=dict(facecolor='black',shrink=0.01,width=0.4,headwidth=3,headlength=3))
    plt.xlabel('Channel width [m]')
    ax_common.set_ylabel('Shelf elevation [m]',labelpad=12)
    
    y_twinx = ax_common.twinx()   
    y_twinx.set_ylabel('Hydrostatic dev. [m]', color= red ,labelpad=23)
    y_twinx.set_yticks([])
    
    
    # Spines for axis
    ax_common.tick_params(labelcolor='None',top=False,bottom=False,left=False,right=False)
    y_twinx.spines['right'].set_visible(False)
    y_twinx.spines['left'].set_visible(False)
    y_twinx.spines['top'].set_visible(False)
    y_twinx.spines['bottom'].set_visible(False)
  
    plt.subplots_adjust(top=0.95,hspace=0)
    
    
  
    
    # Save figures 
    path = str('plots/03_results/01_hydrostaticdeviation')
       
    fname_png = str('hydrostatic_dev_2d_' + str(original_width*2) + '.png')
    fname_pdf = str('hydrostatic_dev_2d_' + str(original_width*2) + '.pdf')
   
    plt.savefig(path + fname_png, format = 'png',dpi=1000,bbox_inches = 'tight')
    plt.savefig(path + fname_pdf, format = 'pdf',dpi=1000,bbox_inches = 'tight')       
        
    plt.show()
