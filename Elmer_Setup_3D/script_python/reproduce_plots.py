#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 17:24:55 2019

@author: jloos
"""

from main import ModelRun

# 2-D
from __plot_hd_2d_line import plot_hd_line_2d
from __Plot_scalar_2d import Plot_scalar_2D
from __plot_profile_velos_2d import plot_profile_velos_2d
from __plot_profile_derivative_2d import plot_profile_derivative_2d
from __plot_abutment_peakdeviation import plot_abutment_peakdeviation
from __plot_rms_2d import plot_rms_2d
from __plot_channelclosing import plot_channelclosing

# 3-D
from __Plot_scalar_3d import Plot_scalar_3d
from __plot_rms_3d import plot_rms_3d

# Appendix
from plot_zs_zb import plot_zs_zb
from plot_zs_derivatives import plot_zs_derivative


def main():
    
    """
    Calls main function to reproduce all plots in a row.   
    """

    
    # Call all plots for figures in the thesis.
    # Uncomment if you do not want to run the function
    
    
    # Figure 5
    plot_hd_line_2d(t1=10,t2=40,width=50000)  
    
    # Figure 6
    Plot_scalar_2D(t1=10,t2=40,width=50000,scalar = "velocity_x")
    Plot_scalar_2D(t1=10,t2=40,width=100000,scalar = "velocity_x")
    
    # Figure 7
    Plot_scalar_2D(t1=10,t2=40,width=50000,scalar = "velocity_y")
    Plot_scalar_2D(t1=10,t2=40,width=100000,scalar = "velocity_y")
    
    # Figure 8
    plot_profile_velos_2d(t1 = 10, width_1 = 50000, width_2 = 100000)
    
    # Figure 9
    Plot_scalar_2D(t1=10,t2=40,width=50000)
    Plot_scalar_2D(t1=10,t2=40,width=50000)
    
    # Figure 10
    plot_profile_derivative_2d(t1 = 40, width = 50000, letter = 'a')
    plot_profile_derivative_2d(t1 = 40, width = 100000, letter = 'b')
    
    # Figure 11
    plot_abutment_peakdeviation(t_end = 100)
    
    # Figure 12
    plot_rms_2d(t1=4,t2=10,t3=40,t4=50)
    
    # Figure 13
    plot_channelclosing(t_end = 100)
    
    # Figure 14
    Plot_scalar_3d(t1=200,amp=250,width=100)
    Plot_scalar_3d(t1=200,amp=250,width=400)
    
    # Figure 15
    plot_rms_3d()
    
    # Figure Append. 18
    plot_zs_zb(t1=100)
    plot_zs_zb(t1=200)
    
    # Figure Append. 19
    plot_zs_derivative(t1=100)
    plot_zs_derivative(t1=200)
    
    
if __name__ == '__main__':
    
    main()