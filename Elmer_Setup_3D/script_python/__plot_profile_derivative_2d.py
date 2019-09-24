# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 15:55:13 2018
"""


from main import ModelRun
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from scipy.interpolate import griddata
import pandas as pd
from __plot_params import params_horizontal
from scipy.interpolate import interp1d
from numpy import diff
from scipy import ndimage
import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

    

def plot_profile_derivative_2d(t1 = None,
                               width = None,
                               prop = None,
                               x_prof = None,
                               letter = None):
                
    '''
    
    Function: plot_profile_derivative_2d
    ---------------
    This function plots for a given timestep "t1", transect "x_prof" and channel
    width "width" a bridging profile and the deviation of hydrostatic equilibrium
    in the first subplot. The second subplots describes the derivative of the
    bridging profile, the deviation of hydrostatic equilibrium and the top 
    surface zs.
    
    Parameters
    ----------
    t1 : int
        timestep for plot
    width : int (w_mod)
        chosen width for plot
    prop : int (0 for regular domain)
        timestep for plot  
    x_prof : int
        transect of the ice-shelf for the 'n'th column (15 layers)
        default is 7, which is approximately -124 m, closely above the channel 
        apex
    letter : str
        annotation of subplots (in general not important and is only used for
        thesis figure)
        
    '''
    
    
    
    # Choose default timesteps if timesteps are not given
    if t1 is None:
            t1 = 50
    # Is set to zero for regular domain size 
    if prop is None:
        prop=0
    
    if x_prof is None:
        x_prof = 7
   
    if letter is None:
        letter = 'a'      
        
    

    # Colors
    orange = '#D55E00'  
    red = '#CC2529'
    wheat = '#948B3D'


                
    #calculation of real width
    original_halfwidth = int(round(np.sqrt(width)*2))
    original_width = original_halfwidth
    


    nrow = 2; ncol = 1;
    fig, axsi = plt.subplots(nrow,ncol,sharex = True, squeeze = False, constrained_layout=True)
    
    ax1 = axsi[0][0]
    ax2 = axsi[1][0]
    
    plt.subplots_adjust(hspace=0.1)
    
    

    
     # Custom params load from __plot_params
    plt.rcParams.update(params_horizontal) 
    

    mr = ModelRun(150,width,0,prop,t1,"2")
    ht = mr.compute_hydrostatic_thickness()

    upper=ht[2]
    lower=ht[3]
    
    scalar_real = mr.get_scalar_real(' sxy')
    
    x_mat = scalar_real[0]
    y_mat = scalar_real[1]
    scalar_mat = scalar_real[2]
    
    
    # compute deviation regular
    upper=ht[2]
    lower=ht[3]
    points=ht[4]
    x_line = ht[0]
    lower = ht[3]
    
    # Calculated ydrostatic thickness
    calc_thickness_bs = ht[1]
    
    # Deviation of hydrostatic equilibrium
    hydrostatic_deviation = calc_thickness_bs - lower
    
    # Plot bridging
    ax1.plot(x_mat[x_prof,:],scalar_mat[x_prof,:],linestyle = '-',linewidth = 1, color = 'black')
    
    # Annotation for cross section
    ax1.text(-2500,0.031,"A",ha = 'left', va='bottom',color = 'y',size=15)
    ax1.text(2500,0.031,"A'",ha = 'right', va='bottom', color = 'y',size=15)
    
    # Vertical lines
    ax1.axvline(0,linestyle='--',linewidth=0.6,color='red')
    ax1.axhline(0,linestyle='--',linewidth=0.6,color='red')
    ax1.tick_params(direction='in',length=3,width=1)
    ax1.set_xticks([])
    ax1.set_xticklabels([])
    ax1.set_xlim(-2500,2500)
    ax1.set_ylim(-0.031,0.031)
    
    ax1.set_ylabel('$\sigma_{xy}$ [$MPa$]', visible = True)
    
    legend_ax = ax1.legend(['Bridging profile $(f(x))$: along x-axis @ y =  ' + str(int(np.mean(y_mat[x_prof,:]))) + ' m'], loc= "upper right")
    frame_ax = legend_ax.get_frame()
    frame_ax.set_facecolor('0.7')
    frame_ax.set_edgecolor('0.7')


    
    # place text box in upper left in axes coords      
    props = dict(boxstyle='round', facecolor='wheat')
    
    
    # Twin plot
    axs = ax1.twinx()
    axs.tick_params(direction='in',length=3,width=1)
    axs.set_xticks([])
    axs.grid(False)
    
    
    axs.set_yticks([-30,-10,0,10,30])
    axs.set_ylim(-40,40)
    axs.set_xticklabels([])
    plt.setp(axs.get_yticklabels(),fontweight = 'bold',color=red)
    axs.plot(x_line, hydrostatic_deviation, linestyle = '-', color = red,linewidth=1)


    title = str('Bridging profile vs. dev. of hydrostatic eq.')
    
    
    axs.text(0.02, 0.95, str('(') + letter + str('.1)') + '  t = ' + str(t1*5) + 'a\n @ cw = ' + str(original_width*2) + 'm', transform=axs.transAxes,
        verticalalignment='top', bbox=props,weight='bold')
    ax2.text(0.02, 0.85, str('(') + letter + str('.2)'), transform=ax2.transAxes,
        verticalalignment='top', bbox=props,weight='bold')
       
    ax2.set_xlabel('Across flow distance [m]') 
    ax2.set_ylabel('$\\frac{dy}{dx}\:[]$') 
    
    
    # Calculate first diravitive
    dydx = diff(scalar_mat[x_prof,:])/diff(x_mat[x_prof,:])
    new_x = x_mat[x_prof,:]
    dx = new_x[1]-new_x[0]
    dxdx = dx**2
    df = np.diff(scalar_mat[x_prof,:]) /dx
    gf = ndimage.gaussian_filter1d(scalar_mat[x_prof,:],sigma = 1, order=1, mode='wrap') /dx
    
    
    # Plot for zs
    ax_upper = fig.add_subplot(212)
    ax_upper.plot(x_line,upper)
    ax_upper.set_axis_off()
    ax_upper.set_xlim(-2500,2500)
    legend_ax_upper = ax_upper.legend(['$z_{s}$'], loc= "upper right")
    frame_ax_upper = legend_ax_upper.get_frame()
    frame_ax_upper.set_facecolor('0.7')
    frame_ax_upper.set_edgecolor('0.7')
    
    
    
    
    #dydx = np.append(dydx,1)
    ax2.plot(x_mat[x_prof,:],-gf,color = 'black',linestyle=':',linewidth = 1) 
    
    
    ax2.set_xticks([-1500,-1000,0,1000,1500])
    ax2.set_xticklabels([-1500,-1000,0,1000,1500])
    ax2.set_xlim(-2500,2500)
    ax2.tick_params(direction='in',length=3,width=1)
    
    legend_ax3 = ax2.legend(['$-\\frac{d}{dx} f(x)$'],loc="lower right", prop=dict(weight='bold'))
    frame_ax3 = legend_ax3.get_frame()
    frame_ax3.set_facecolor('0.7')
    frame_ax3.set_edgecolor('0.7')

    # Plot hydrostatic deviation
    ax3 = ax2.twinx()
    ax3.plot(x_line, hydrostatic_deviation, linestyle= '-', color = red,linewidth=1)
    plt.setp(ax3.get_yticklabels(),fontweight = 'bold',color=red)
    ax3.set_yticks([-10,-5,0,5,10])
    ax3.tick_params(direction='in',length=3,width=1)
    ax3.grid(False)
    
    fig.text(1,0.5,'Hydrostatic dev. [m]',va = 'center',rotation = 'vertical',color=red)

   
    # Save figures    
    path = str('plots/03_results/03_shearstress/')
       
    fname_png = str('profile_deviation_dev_2d_' + str(original_width*2) + '.png')
    fname_pdf = str('profile_deviation_dev_2d_' + str(original_width*2) + '.pdf')

    
    plt.savefig(path + fname_png, format = 'png',dpi=1000,bbox_inches='tight')
    plt.savefig(path + fname_pdf, format = 'pdf',dpi=1000,bbox_inches='tight')       
        
    plt.show()    
    
    
