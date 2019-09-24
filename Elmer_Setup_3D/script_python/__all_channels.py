from main import ModelRun
from __plot_params import params_bridging_2d,params_horizontal


from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
import matplotlib.font_manager
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy.interpolate import griddata,interpolate

 
# List of widths
list_widths_1 = list(np.arange(10000,100000,10000))
list_widths_2 = list(np.arange(200000,1000000,100000))
          
            
list_widths_1.extend(list_widths_2)

# Define rc_params for figure
plt.rcParams.update(params_horizontal)   

fig, ax1 = plt.subplots()  
for width in list_widths_1:
       
    
    # Empty lists for plot
    hd_list_tight = []
    hd_list_wide = []
    hd_list_extrawide = []
    
    sxy_list = []
    sxy_list_extent = []
    sxy_list_wideextent = []
    
    channel_abuts_1 = []
    channel_abuts_2 = []
    channel_abuts_3 = []
    
    ht = ModelRun(150,width,0,0,0,"2").compute_hydrostatic_thickness() 
    upper=ht[2]
    lower=ht[3]
    new_x = ht[0]
    ax1.set_xlim(-2500,2500)
    
    ax1.plot(new_x,lower,'k-',linewidth=0.5) 
    fig.add_subplot(ax1)    
    ax1.set_xlabel('Channel width [m]')
    ax1.set_ylabel('Channel height [m]')
    
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    
    ax1.tick_params(direction='in',length=4,width=2)
    ax1.set_yticks([-260,-200,-160])
    ax1.set_xticks([-2000,-500,0,500,2000])
    plt.setp(ax1.get_xticklabels(),fontweight = 'bold')
    plt.setp(ax1.get_yticklabels(),fontweight = 'bold')



# save figures
plt.savefig('plots/02_methods/all_channels.png', format = 'png',dpi=1000) 
plt.savefig('plots/02_methods/all_channels.pdf', format = 'pdf',dpi=1000)   
plt.show()  