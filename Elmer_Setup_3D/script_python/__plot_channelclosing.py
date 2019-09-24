# Custom modules
from main import ModelRun
from __plot_params import params_horizontal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'


def plot_channelclosing(t_end=None):   
    
    '''
      
    Function: plot_channelclosing
    -----------------------------
    This functions generates a figure containing six subplots with different
    widths (list_widths) for a given time range. Each subplot depicts the 
    original channel amplitude elevation.
    
    Parameters
    ----------
    t_end : int
        timestep        
        
    '''
   
   
    # Get domain-size
    points = ModelRun(150,20000,0,0,50,"2").Points
    points_extent =  ModelRun(150,20000,0,0,50,"2").Points
    points_wideextent =  ModelRun(150,20000,0,0,50,"2").Points
    
    x = points[:,0]
    x_extent = points_extent[:,0]
    x_wideextent = points_wideextent[:,0]
    
     # List of widths for extended and tight domain
    list_widths = list(np.arange(20000,100000,10000))
    list_widths = [50000,80000,90000,300000,500000,700000]
    
    # Numeration for subplots
    nums = ["(a)","(b)","(c)","(d)","(e)","(f)"]
    

    fig, axs = plt.subplots(2,3, sharex=True, sharey = True) 
    
    
    # Custom model load from __plot_params
    plt.rcParams.update(params_horizontal)
    plt.subplots_adjust(hspace=0,wspace=0)
    
    for width,ax1,num in zip(list_widths,axs.flat,nums):
       
               
        
        # Empty lists for plot        
        channel_abuts_1 = []
        channel_abuts_2 = []
        channel_abuts_3 = []
        
        time = []
        
        
        if t_end is None:
            t_end=200
        for i in range(3,t_end):
            hd_regular = ModelRun(150,width,0,0,i,"2").compute_hydrostatic_thickness()
            hd_extended = ModelRun(150,width,0,'extent',i,"2").compute_hydrostatic_thickness()
            hd_extended_wide = ModelRun(150,width,0,'wideextent',i,"2").compute_hydrostatic_thickness()
                       
            channel_abuts_1.append(hd_regular[10])
            channel_abuts_2.append(hd_extended[10])
            channel_abuts_3.append(hd_extended_wide[10])
            
            
            time.append(i*5)
            
          
            
        #  Make plot, cut frame, properties defined in rc.params
        def original_width(width):
            original_width = int(round(np.sqrt(width)*2)*2)
            return original_width                   
        
        
        
        # place text box in upper left in axes coords      
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)        
        ax1.text(0.05, 0.92, num + '  width = ' + str(original_width(width)) + 'm', transform=ax1.transAxes, 
                verticalalignment='top', bbox=props, weight='bold')
        
        

        ax1.tick_params(direction='in',length=3,width=1)

        ax1.get_xaxis().set_ticks([50, 250, 450])
        ax1.yaxis.set_major_locator(MaxNLocator(3))

        ax1.set_xticks([50,250,450])
      
        ax1.plot(time, channel_abuts_1, linestyle='-', color='black',linewidth = 1)
        ax1.plot(time, channel_abuts_2, linestyle='--', color='black', linewidth = 1)
        ax1.plot(time, channel_abuts_3, linestyle=':',color='black', linewidth = 1)
        ax1.tick_params(direction='in',length=3,width=1)
        
        fig.add_subplot(ax1)
        
        
    fig.text(0,0.5,'Elevation of channel amp. [m]',va = 'center',rotation = 'vertical',fontsize=6.5)
    fig.text(0.49,0.05,'Time [a]',va = 'center',fontsize=6.5)

    legend = ax1.legend(['regular', 'extended', 'wide-extended'],loc='center right')                
    frame = legend.get_frame()
    frame.set_facecolor('0.7')
    frame.set_edgecolor('0.7')
    
    # Save figures
    path = str('plots/03_results/06_channelclosing/')
    fname_eps = str('channelclosing_2d_' + str(t_end*5) + 'a' + '.pdf')
    fname_png = str('channelclosing_2d_' + str(t_end*5) + 'a' + '.png')
    
    fig.savefig(path + fname_eps, format = 'pdf',dpi=1000,bbox_inches = 'tight') 
    fig.savefig(path + fname_png, format = 'png',dpi=1000,bbox_inches = 'tight') 
    
    
    plt.show()