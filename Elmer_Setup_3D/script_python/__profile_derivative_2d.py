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

    

def profile_derivative_2d(t1 = None,
                          width = None,
                          prop = None,
                          x_prof= None):
                

    
    
    
        # Choose default timesteps if timesteps are not given
    if t1 is None:
            t1 = 50

     
    if prop is None:
        prop=0
    
    if x_prof is None:
        x_prof = 7
         
    



    # Colors
    orange = '#D55E00'  
    red = '#CC2529'
    wheat = '#948B3D'


                
   
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
    calc_thickness_bs = ht[1]
    

    
    hydrostatic_deviation = calc_thickness_bs - lower
    
    ax1.plot(lower)
    ax1.plot(x_mat[x_prof,:],scalar_mat[x_prof,:],linestyle = '-',linewidth = 1, color = 'black')
    

    
    ax1.axvline(0,linestyle='--',linewidth=0.6,color='red')
    ax1.axhline(0,linestyle='--',linewidth=0.6,color='red')
    ax1.tick_params(direction='in',length=3,width=1)
    ax1.set_xticks([])
    ax1.set_xticklabels([])
    ax1.set_xlim(-2500,2500)
    ax1.set_ylim(-0.031,0.031)
    
    ax1.set_ylabel('$\sigma_{xy}$ [$MPa$]', visible = True)
    
    legend_ax = ax1.legend(['Bridging profile $f(x)$: along x-axis @ y =  ' + str(int(np.mean(y_mat[x_prof,:]))) + ' m'], loc= "upper right")
    frame_ax = legend_ax.get_frame()
    frame_ax.set_facecolor('0.7')
    frame_ax.set_edgecolor('0.7')


    
    # place text box in upper left in axes coords      
    props = dict(boxstyle='round', facecolor='wheat')
    
    
  
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
    
    
    axs.text(0.02, 0.95, str('(b.1)') + '  t = ' + str(t1*5) + 'a\n @ cw = ' + str(original_width*2) + 'm', transform=axs.transAxes,
        verticalalignment='top', bbox=props,weight='bold')
    ax2.text(0.02, 0.95, str('(b.2)'), transform=ax2.transAxes,
        verticalalignment='top', bbox=props,weight='bold')
       
    ax2.set_xlabel('Across flow distance [m]') 
    ax2.set_ylabel('$\\frac{dy}{dx}$') 
    
    
    # Calculate first diravitive
    dydx = diff(scalar_mat[x_prof,:])/diff(x_mat[x_prof,:])
    new_x = x_mat[x_prof,:]
    dx = new_x[1]-new_x[0]
    dxdx = dx**2
    df = np.diff(scalar_mat[x_prof,:]) /dx
    gf = ndimage.gaussian_filter1d(scalar_mat[x_prof,:],sigma = 1, order=1,mode='wrap') /dx
    
    
    
    
    
    #dydx = np.append(dydx,1)
    ax2.plot(x_mat[x_prof,:],-gf,color = 'black',linestyle=':',linewidth = 1)    
    
    
    ax2.set_xticks([-1500,-1000,0,1000,1500])
    ax2.set_xticklabels([-1500,-1000,0,1000,1500])
    #ax2.set_ylim(-0.00025,0.00025)
    ax2.set_xlim(-2500,2500)
    #ax2.set_xticklabels([])
    ax2.tick_params(direction='in',length=3,width=1)
    
   
    legend_ax3 = ax2.legend(['$\\frac{d}{dx} f(x)$'],loc="upper right", prop=dict(weight='bold'))
    frame_ax3 = legend_ax3.get_frame()
    frame_ax3.set_facecolor('0.7')
    frame_ax3.set_edgecolor('0.7')
    
   
    
    ax3 = ax2.twinx()
    ax3.plot(x_line, hydrostatic_deviation, linestyle= '-', color = red,linewidth=1)
    #ax3.set_ylim(-40,40)

    plt.setp(ax3.get_yticklabels(),fontweight = 'bold',color=red)
    ax3.set_yticks([-10,0,5])
    #fig.suptitle(title + ' @ cw = ' + str(original_width*2) + 'm', weight='bold', fontsize = 7, y = 0.95)
    ax3.tick_params(direction='in',length=3,width=1)
    ax3.grid(False)
   
    fig.text(1,0.5,'Hydrostatic dev. [m]',va = 'center',rotation = 'vertical',color=red)

   
    # Save figures    
    path = str('plots/Final_plots/')
       
    fname_png = str('TESTprofile_deviation_dev_2d_' + str(original_width*2) + '.png')
    fname_pdf = str('TESTprofile_deviation_dev_2d_' + str(original_width*2) + '.pdf')

    
    plt.savefig(path + fname_png, format = 'png',dpi=1000,bbox_inches='tight')
    plt.savefig(path + fname_pdf, format = 'pdf',dpi=1000,bbox_inches='tight')       
        
    plt.show()    
    
    
def profile_multiple_derivative_2d(t1 = None,
                          prop = None,
                          x_prof= None):
                

    
    
    
        # Choose default timesteps if timesteps are not given
    if t1 is None:
        t1 = 50
    
    if prop is None:
        prop=0
    
    if x_prof is None:
        x_prof = 7
         
    



    # Colors
    orange = '#D55E00'  
    red = '#CC2529'
    wheat = '#948B3D'

    list_widths_1 = list(np.arange(10000,100000,10000))
    list_widths_2 = list(np.arange(200000,900000,100000))
    
    list_widths_1.extend(list_widths_2)
                
   
    
    


    nrow = 1; ncol = 1;
    fig, axsi = plt.subplots(nrow,ncol,sharex = True, squeeze = False, constrained_layout=True)
    
    
    
    ax2 = axsi[0][0]
    

    
     # Custom params load from __plot_params
    plt.rcParams.update(params_horizontal) 
    
    for width in list_widths_2:
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
        calc_thickness_bs = ht[1]
        
        
        
        hydrostatic_deviation = calc_thickness_bs - lower
        original_halfwidth = int(round(np.sqrt(width)*2))
        original_width = original_halfwidth
        
        
        
          
        
        
        
        title = str('Bridging profile vs. dev. of hydrostatic eq.')
        
        
         # place text box in upper left in axes coords      
        props = dict(boxstyle='round', facecolor='wheat')
        ax2.text(0.02, 0.95, str('(a)'), transform=ax2.transAxes,
            verticalalignment='top', bbox=props,weight='bold')
           
        ax2.set_xlabel('Across flow distance [m]') 
        ax2.set_ylabel('$\\frac{dy}{dx}$') 
        
        
        # Calculate first diravitive
        dydx = diff(scalar_mat[x_prof,:])/diff(x_mat[x_prof,:])
        new_x = x_mat[x_prof,:]
        dx = new_x[1]-new_x[0]
        dxdx = dx**2
        df = np.diff(scalar_mat[x_prof,:]) /dx
        gf = ndimage.gaussian_filter1d(scalar_mat[x_prof,:],sigma = 1, order=1,mode='wrap') /dx
        
        
        
        
        
        #dydx = np.append(dydx,1)
        ax2.plot(x_mat[x_prof,:],-gf,linestyle='-',linewidth = 0.5)    
        
           
        legend_ax3 = ax2.legend(['$\\frac{d}{dx} f(x)$'],loc="upper right", prop=dict(weight='bold'))
        frame_ax3 = legend_ax3.get_frame()
        frame_ax3.set_facecolor('0.7')
        frame_ax3.set_edgecolor('0.7')
        
           
        
        ax3 = ax2.twinx()
         ax3.set_yticklabels([])
        ax3.plot(x_line, hydrostatic_deviation, linestyle= '-',linewidth=0.5)
#        
        

    plt.setp(ax3.get_yticklabels(),fontweight = 'bold',color=red)
    ax3.set_yticklabels([])
    ax3.set_yticks([-10,0,5])
    #fig.suptitle(title + ' @ cw = ' + str(original_width*2) + 'm', weight='bold', fontsize = 7, y = 0.95)
    ax3.tick_params(direction='in',length=3,width=1)
    ax3.grid(False)
    ax2.set_xticks([-1500,-1000,0,1000,1500])
    ax2.set_xticklabels([-1500,-1000,0,1000,1500])
    #ax2.set_ylim(-0.00025,0.00025)
    ax2.set_xlim(-2500,2500)
    #ax2.set_xticklabels([])
    ax2.tick_params(direction='in',length=3,width=1)
    fig.text(1,0.5,'Hydrostatic dev. [m]',va = 'center',rotation = 'vertical',color=red)

   
    # Save figures    
    path = str('plots/Final_plots/')
       
    fname_png = str('multiple_profile_deviation_dev_2d_' + str(original_width*2) + '.png')
    fname_pdf = str('multiple_profile_deviation_dev_2d_' + str(original_width*2) + '.pdf')

    
    plt.savefig(path + fname_png, format = 'png',dpi=1000,bbox_inches='tight')
    plt.savefig(path + fname_pdf, format = 'pdf',dpi=1000,bbox_inches='tight')       
        
    plt.show()    