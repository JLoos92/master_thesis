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
from __plot_params import params_horizontal,params_vertical,params_horizontal_single


from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
import matplotlib.font_manager
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy.interpolate import griddata,interpolate




def scalar_grid_plot(t1 = None,
                     t2 = None,
                     width = None,
                     prop = None,
                     scalar = None,
                     dim = None):
                
    
    # Choose default timesteps if timesteps are not given
    if (t1,t2) is None:
            t1 = 50
            t2 = 100

    if dim is not None:
        dim = str('j2')         
        print(dim)
    if prop is None:
        prop = 0        
        
    # Time + numbers for ax
    nums = ["(a)","(b)","(c)","(d)"]
    times = [t1, t2]
    
    
    
    # Colors
    orange = '#D55E00'
    blue = '#396AB1'    
    red = '#CC2529'
    wheat = '#948B3D'
    
    # Define plot layout
    nrow = 2; ncol = 1;
    fig, ax1 = plt.subplots(nrows=nrow, ncols=ncol,sharey = True, squeeze = False, constrained_layout = True)
    
    # For common axis add subplot
    ax_common = fig.add_subplot(111,frameon=False) 
    plt.subplots_adjust(top=2)

    plt.tick_params(labelcolor=False,top=False,bottom=False,left=False,right=False)
    
    # Custom params load from __plot_params
    plt.rcParams.update(params_horizontal_single) 
    
    
    for ax,num,t in zip(ax1.reshape(-1),nums,times): 
        
            #Extract data from modelrun
            mr = ModelRun(150,width,0,prop,t,"2")
            ht = mr.compute_hydrostatic_thickness()
            

            # original channel_width
            original_width = int(round(np.sqrt(width)*2)) 
            
            # Choose scalar from run (default is bridging)
            if scalar is None:
                scalar_real = mr.get_scalar_real(' sxy')
                x_mat = scalar_real[0]
                y_mat = scalar_real[1]
                scalar_mat = scalar_real[2]
                
                cmap = str("bwr")
                title = str("Bridging $\sigma_{xy}$ [$MPa$]")
                scalarname = str('bridging')
                vmin = -0.03
                vmax = 0.03
                pass
            
            elif scalar is not None:                
                scalar_real = mr.get_scalar_real(str(scalar))
                x_mat = scalar_real[0]
                y_mat = scalar_real[1]
                scalar_mat = scalar_real[2]
                vmax = 0.3
                vmin = -0.3
                cmap = str("RdBu")
                scalarname = str('velo_x')
                title = str("Velocity $u_{x}$ [$m\: a^{-1}$]")
                
            elif dim is str('y'):     
                cmap = str("RdYlBu")                
                
                scalar_real = mr.get_scalar_real(str(scalar))
                x_mat = scalar_real[0]
                y_mat = scalar_real[1]
                scalar_mat = scalar_real[3]
                vmax = -0.8
                vmin = -0.3
                scalarname = str('velo_y')
                title = str("Velocity $u_{y}$ [$m\: a^{-1}$]")
                
                
     #   elif dim is str('j2'):     
            cmap = str("RdYlBu")                
            
            scalar_real_sxx = mr.get_scalar_real('sxx')
            scalar_real_syy = mr.get_scalar_real(' syy')
            
            x_mat = scalar_real_sxx[0]
            y_mat = scalar_real_sxx[1]
            
            sxx_mat = scalar_real_sxx[2] 
            syy_mat = scalar_real_syy[2]
            
            scalar_mat = 1/2*((sxx_mat**2)+(syy_mat**2))
            
            vmax = 2
            vmin = -2
            scalarname = str('Second invariant')
            title = str("Effective stress [$MPa$]")
                
                
            # Get boundaries
            upper=ht[2]
            lower=ht[3]
            new_x = ht[0]
        
            # Calculated hydrostatic thickness
            calc_thickness_bs = ht[1]
              
            points = ht[4]
            x = points[:,0]
            y = points[:,1]
            
            y = np.flipud(y)

           
            # Points for triangulation
            points = [x,y]
            points = np.asarray(points)
            points = points.transpose()
            
            
            im = ax.pcolormesh(x_mat,y_mat,scalar_mat,cmap = cmap ,vmin = vmin, vmax = vmax, shading = "gouraud")
            
            fig.subplots_adjust(top=1.1)
            
            ax_cbar = fig.add_axes([0.11, 0.86, 0.8, 0.025])
            plt.colorbar(im, cax=ax_cbar, orientation='horizontal',ticklocation = 'top',extend = 'both')
            ax_cbar.xaxis.set_tick_params(pad=1,size=1)
           
            

            # Set labels
            ax.set_xlabel('Channel width [m]', visible = True)
           

            ax.set_xticks([-original_width,0,original_width])
            ax.set_yticks([-250,-200,-100,-50])
            
            plt.setp(ax.get_xticklabels(),fontweight = 'bold')
            plt.setp(ax.get_yticklabels(),fontweight = 'bold')
            
            upper=ht[2]
            lower=ht[3]
            points=ht[4]
            x_line = ht[0]
            
            # Plot lines
            ax.plot(x_line,lower,'k-',linewidth=2)
            ax.fill_between(x_line, lower, lower.min(), color='w')
            #ax.plot(x_line,upper,'b-',linewidth=1)
            ax.fill_between(x_line, upper, upper.max(), color='w')
            ax.tick_params(direction='in',length=4,width=2)
            ax.set_xlim(-1500,1500)
            ax.set_ylim(-290,0)
             
            
            # Hydrostatic thickness
            ax.plot(new_x,calc_thickness_bs,linestyle = '--',color = blue,linewidth=2)
            
            # Legend
            legend_ax2 = ax.legend(['Modelled thickness','Hydrostatic thickness'],loc="upper right", prop=dict(weight='bold'))
            frame_ax2 = legend_ax2.get_frame()
            frame_ax2.set_facecolor('0.7')
            frame_ax2.set_edgecolor('0.7')
            
            
            
            # compute deviation regular
            lower = ht[3]
            calc_thickness_bs = ht[1]


            
            hydrostatic_deviation = calc_thickness_bs - lower
            
            # place text box in upper left in axes coords      
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
            

            
            
            # Second axis
            axs = ax.twinx()
            axs.tick_params(direction='in',length=4,width=2)
              
            axs.set_yticks([-30,-10,0,10,30])
            axs.set_ylim(-40,40)
            
            plt.setp(axs.get_yticklabels(),fontweight = 'bold',color=red)
            
            axs.plot(x_line, hydrostatic_deviation, linestyle = ':',color = red ,linewidth=3)

            ax.text(0.04, 0.92, num + '  t = ' + str(t*5) + 'a', transform=ax.transAxes, 
                verticalalignment='top', bbox=props, weight='bold',fontsize=9)
            
            # Spines in subplot
            axs.spines['top'].set_visible(False)
            #axs.spines['bottom'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            fig.add_subplot(ax)

            fig.add_subplot(ax)
     
        
    fig.suptitle(title + ' @ cw = ' + str(original_width*2) + 'm')


    # Set legend downer off
    ax1[1][0].get_legend().remove()         
    
    
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
  
    plt.subplots_adjust(top=0.85,hspace=0)
        
                 
            





       
    # Save figure   
    path = str('plots/Final_plots/scalar_2d/')
    fname_png = str('j2_dev_2d__' + str(original_width*2) + '.png')
    fname_pdf = str('j2_dev_2d__' + str(original_width*2) + '.pdf')

    plt.savefig(path + str(scalarname) + fname_png, format = 'png',dpi=1000,bbox_inches='tight')
    plt.savefig(path + str(scalarname) + fname_pdf, format = 'pdf',dpi=1000,bbox_inches='tight')
    
    plt.show() 
        
  