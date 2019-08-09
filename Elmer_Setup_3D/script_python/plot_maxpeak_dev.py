# Custom modules
from main import ModelRun
from __plot_params import params_horizontal
from local_minimum import local_minimum_1d
#
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter



def compute_maxpeak_dev(t_end=None,
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
    list_widths = [50000,70000,80000,100000,200000,300000]
    
    width_1 = 50000
    width_2 = 50000  
    width_3 = 50000
    
    # Time + numbers for ax
    nums = ["(a)","(b)","(c)","(d)","(e)","(f)"]
    
    orange = '#D55E00'
    fig, axs = plt.subplots(2,3, sharex=True, sharey = True) 
    
    
     # Custom model load from __plot_params
    plt.rcParams.update(params_horizontal)
    plt.subplots_adjust(hspace=0,wspace=0)
    
    for width,ax1,num in zip(list_widths,axs.flat,nums):
        # Define rc_params for figure
               
        
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
        
        time = []
        
        
        if t_end is None:
            t_end=200
        for i in range(3,t_end):
            hd_regular = ModelRun(150,width,0,0,i,"2").compute_hydrostatic_thickness()
            hd_extended = ModelRun(150,width,0,'extent',i,"2").compute_hydrostatic_thickness()
            hd_extended_wide = ModelRun(150,width,0,'wideextent',i,"2").compute_hydrostatic_thickness()
#           
#            # Get bridging stresses
#            s_xy =  ModelRun(150,width_1,0,0,i,"2").sxy
#            s_xy_extent = ModelRun(150,width_2,0,0,i,"2").sxy
#            s_xy_wideextent = ModelRun(150,width_3,0,0,i,"2").sxy
#            
#            # Regular
#            df = pd.DataFrame({'sxy':s_xy,'x':x})
#            df = df[df.x>0]
#            sxy_array = df.iloc[:,0].values
#            rms_sxy = np.sqrt(np.mean(sxy_array**2))
#            
#            # Extended
#            df_extent = pd.DataFrame({'sxy_extent':s_xy_extent,'x_extent':x_extent})
#            df_extent = df_extent[df_extent.x_extent>0]
#            df_extent = df_extent[df_extent.x_extent<10000]
#            sxy_array_extent = df_extent.iloc[:,0].values
#            rms_sxy_extent = np.sqrt(np.mean(sxy_array_extent**2))
#            #import pdb;pdb.set_trace()
#            
#            # Wide extended
#            df_wideextent = pd.DataFrame({'sxy_wideextent':s_xy_wideextent,'x_wideextent':x_wideextent})
#            df_wideextent = df_wideextent[df_wideextent.x_wideextent>0]
#            df_wideextent = df_wideextent[df_wideextent.x_wideextent<10000]
#            sxy_array_wideextent = df_wideextent.iloc[:,0].values
#            rms_sxy_wideextent = np.sqrt(np.mean(sxy_array_wideextent**2))
#            
#            # Get values of max peak deviation (channel's amplitude) in percent
#            ht_regular = hd_regular[7]
#            ht_extended = hd_extended[7]
#            ht_extended_wide = hd_extended_wide[7]
#            
#            # Append values
#            hd_list_tight.append(ht_regular)
#            hd_list_wide.append(ht_extended)
#            hd_list_extrawide.append(ht_extended_wide)
#            
#            
#            sxy_list.append(rms_sxy)
#            sxy_list_extent.append(rms_sxy_extent)
#            sxy_list_wideextent.append(rms_sxy_wideextent)
            
            
#            peak_1 = local_minimum_1d(width_1,t)
#            peak_2 = local_minimum_1d(width_2,t)
#            peak_3 = local_minimum_1d(width_2,t)
            
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
        
        ax1.text(0.04, 0.95, num + '  width = ' + str(original_width(width)) + 'm', transform=ax1.transAxes, 
                verticalalignment='top', bbox=props, weight='bold',fontsize=5)
        
        
        
#        #ax1.set_frame_on(frameon = False)
#        ax1.plot(time,hd_list_tight,'k-',linewidth= 1)
#        ax1.plot(time,hd_list_wide,'k--',linewidth= 1)
#        ax1.plot(time,hd_list_extrawide,'k:',linewidth= 1)
    
        
        
#        
#        #ax1.title('The maximum peak deviation', y = 1.05)
#        ax1.set_xlim(15,t_end*5)
#        ax1.set_ylim(0,25)       
#        #ax1.set_xlabel('Time [a]',labelpad=5)
        #ax1.set_ylabel('Max. peak deviation [%]',labelpad=5)
        ax1.tick_params(direction='in',length=3,width=1.5)
        
#        
#        ax1.get_yaxis().set_ticks([5,10,15,20])
        ax1.get_xaxis().set_ticks([50, 250, 450])
        ax1.yaxis.set_major_locator(MaxNLocator(3))
        
        
        #plt.setp(ax1.get_xticklabels(),fontweight = 'bold')
        #plt.setp(ax1.get_yticklabels(),fontweight = 'bold')
        ax1.set_xticks([50,250,450])
#        ax1.set_yticks([5,10,15,20])
        #ax1.set_yticklabels([5,15,20])
#        # Plt properties for ax2 (bridging)  
#        ax2 = ax1.twinx()
#        
        ax1.plot(time, channel_abuts_1, linestyle='-',color=orange)
        ax1.plot(time, channel_abuts_2, linestyle='--',color=orange)
        ax1.plot(time, channel_abuts_3, linestyle=':',color=orange)
        #ax2.set_ylabel('Abs. elevation of channel [m]',labelpad = 15, color = 'b')
        ax1.tick_params(direction='in',length=3,width=1)
        
#  
#        # Fontweight ax2
#        plt.setp(ax1.get_xticklabels(),fontweight = 'bold')
#        plt.setp(ax1.get_yticklabels(),fontweight = 'bold')
        
        #fig.text(0.98,0.5,'Elevation of channel amp. [m]',va = 'center',rotation = 'vertical',color=orange,fontsize=7)
        fig.text(0,0.5,'Elevation of channel amp. [m]',va = 'center',rotation = 'vertical')
        fig.text(0.49,0,'Time [a]',va = 'center')
        
        fig.add_subplot(ax1)
        
        


    legend = ax1.legend(['regular', 'extended', 'wide-extended'],loc='center right')
                
    frame = legend.get_frame()
    frame.set_facecolor('0.7')
    frame.set_edgecolor('0.7')
#    
    path = str('plots/Final_plots/max_peak_dev_2d/')
    fname_eps = str('maxmodel_2d_' + str(width_1)+ str(width_2) + '_' + str(t_end*5) + 'a' + '.pdf')
    fname_png = str('maxmodel_2d'+ str(width_1)+ str(width_2) + '_' + str(t_end*5) + 'a' + '.png')
    
    fig.savefig(path + fname_eps, format = 'pdf',dpi=1000,bbox_inches = 'tight') 
    fig.savefig(path + fname_png, format = 'png',dpi=1000,bbox_inches = 'tight') 
    
    
    