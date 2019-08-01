mr_wideextent = ModelRun(150,20000,0,'wideextent',50,'2')
ht_wideextent = mr_wideextent.compute_hydrostatic_thickness()

import pandas as pd

# Custom modules
from main import ModelRun
from __plot_params import params

#
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter




def compute_maxpeak_dev(t=None,
                        width=None,
                        **kwargs):   
    
    '''
      
    Function: compute_maxpeak_dev 
    
    Parameters
    ----------
    t : int
        timestep        
    
    width : kwarg
        width of a channel
      
    
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
    list_widths = [50000]
    
    width_1 = 50000
    width_2 = 100000  
    width_3 = 500000
    
    for width in list_widths:
        # Define rc_params for figure
        fig, ax1 = plt.subplots()        
        
        # Empty lists for plot
        hd_list_tight = []
        hd_list_wide = []
        hd_list_extrawide = []
        
        sxy_list = []
        sxy_list_extent = []
        sxy_list_wideextent = []
        
        time = []
        
        
        if t is None:
            t=200
        for i in range(3,t):
            hd_regular = ModelRun(150,width_1,0,0,i,"2").compute_hydrostatic_thickness()
            hd_extended = ModelRun(150,width_2,0,0,i,"2").compute_hydrostatic_thickness()
            hd_extended_wide = ModelRun(150,width_3,0,0,i,"2").compute_hydrostatic_thickness()
           
            # Get bridging stresses
            s_xy =  ModelRun(150,width_1,0,0,i,"2").sxy
            s_xy_extent = ModelRun(150,width_2,0,0,i,"2").sxy
            s_xy_wideextent = ModelRun(150,width_3,0,0,i,"2").sxy
            
            # Regular
            df = pd.DataFrame({'sxy':s_xy,'x':x})
            df = df[df.x>0]
            sxy_array = df.iloc[:,0].values
            rms_sxy = np.sqrt(np.mean(sxy_array**2))
            
            # Extended
            df_extent = pd.DataFrame({'sxy_extent':s_xy_extent,'x_extent':x_extent})
            df_extent = df_extent[df_extent.x_extent>0]
            df_extent = df_extent[df_extent.x_extent<10000]
            sxy_array_extent = df_extent.iloc[:,0].values
            rms_sxy_extent = np.sqrt(np.mean(sxy_array_extent**2))
            #import pdb;pdb.set_trace()
            
            # Wide extended
            df_wideextent = pd.DataFrame({'sxy_wideextent':s_xy_wideextent,'x_wideextent':x_wideextent})
            df_wideextent = df_wideextent[df_wideextent.x_wideextent>0]
            df_wideextent = df_wideextent[df_wideextent.x_wideextent<10000]
            sxy_array_wideextent = df_wideextent.iloc[:,0].values
            rms_sxy_wideextent = np.sqrt(np.mean(sxy_array_wideextent**2))
            
            # Get values of max peak deviation (channel's amplitude) in percent
            ht_regular = hd_regular[7]
            ht_extended = hd_extended[7]
            ht_extended_wide = hd_extended_wide[7]
            
            # Append values
            hd_list_tight.append(ht_regular)
            hd_list_wide.append(ht_extended)
            hd_list_extrawide.append(ht_extended_wide)
            
            
            sxy_list.append(rms_sxy)
            sxy_list_extent.append(rms_sxy_extent)
            sxy_list_wideextent.append(rms_sxy_wideextent)
                    
            time.append(i*5)
            
          
            
        #  Make plot, cut frame, properties defined in rc.params
        def original_width(width):
            original_width = int(round(np.sqrt(width)*2)*2)
            return original_width                   
        
        #ax1.set_frame_on(frameon = False)
        ax1.plot(time,hd_list_tight,'k-',linewidth= 2)
        ax1.plot(time,hd_list_wide,'k--',linewidth= 2)
        ax1.plot(time,hd_list_extrawide,'k:',linewidth= 2)
        legend = ax1.legend(['cw = ' + str(original_width(width_1))+' m', 'cw = ' + str(original_width(width_2))+' m', 'cw = ' + str(original_width(width_3))+' m)'],loc=1)
        #legend = ax1.legend(['cw = ' + str(original_width)+' m (regular domain)'],loc=2)        
        
                
        frame = legend.get_frame()
        frame.set_facecolor('0.7')
        frame.set_edgecolor('0.7')
        
        
        # Custom model load from __plot_params
        plt.rcParams.update(params) 
        #ax1.title('The maximum peak deviation', y = 1.05)
        ax1.set_xlim(15,t*5)
        ax1.set_ylim(0,25)       
        ax1.set_xlabel('Time [a]',labelpad=5)
        ax1.set_ylabel('Max. peak deviation [%]',labelpad=5)
        ax1.tick_params(direction='in',length=4,width=2)
        
        
        ax1.get_yaxis().set_ticks([5,10,15,20])
        ax1.get_xaxis().set_ticks([100,200,300,400])
        
        plt.setp(ax1.get_xticklabels(),fontweight = 'bold')
        plt.setp(ax1.get_yticklabels(),fontweight = 'bold')
        
#        # Plt properties for ax2 (bridging)  
#        ax2 = ax1.twinx()
#        
#        ax2.plot(time, sxy_list, 'b-')
#        #ax2.plot(time, sxy_list_extent, 'b--')
#        #ax2.plot(time, sxy_list_wideextent, 'b:')
#        ax2.set_ylabel('$\sigma_{xy}$ [MPa]',labelpad = 15, color = 'b')
#        
#        for tl in ax2.get_yticklabels():
#            tl.set_color('b')
#            
#        ax2.set_ylim(0.0010,0.0045)  
#        ax2.tick_params('y')
        
#        
#        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
#        
#        ax1.yaxis.set_major_locator(MaxNLocator(prune='both'))
#        ax2.yaxis.set_major_locator(MaxNLocator(prune='both'))
#        
        
        #legend_ax2 = ax2.legend(['Bridging (regular domain)','Bridging (extended domain)','Bridging (wide extended domain)'],loc=1)
#        legend_ax2 = ax2.legend(['Bridging (regular domain)'],loc=1)
#        frame_ax2 = legend_ax2.get_frame()
#        frame_ax2.set_facecolor('0.7')
#        frame_ax2.set_edgecolor('0.7')
#        
        #plt.title('Peak deviation and bridging stresses for a channel width of ' + str(original_width) + ' m', y = 1.05)
        
        path = str('plots/')
        fname_eps = str('maxpeak_dev_2d__all' + str(original_width) + '_' + str(t*5) + 'a' + '.eps')
        fname_png = str('maxpeak_dev_2d__all' + str(original_width) + '_' + str(t*5) + 'a' + '.png')
        
        fig.savefig(path + fname_eps, format = 'eps',dpi=1000,bbox_inches = 'tight') 
        fig.savefig(path + fname_png, format = 'png',dpi=1000,bbox_inches = 'tight') 