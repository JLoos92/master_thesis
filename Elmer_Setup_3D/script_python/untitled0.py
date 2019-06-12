#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:33:32 2019

@author: jloos
"""

from __main__ import ModelRun



hd_list_tight = []
hd_list_wide = []
hd_list_ultra_wide = []

# Deviation of max. peaks of the amp.

for i in range(3,200):
    hd_tight = ModelRun(150,23000,0,0,i,"2").compute_hydrostatic_thickness()
    hd_wide = ModelRun(150,2200,0,0,i,"2").compute_hydrostatic_thickness()
    hd_ultra_wide = ModelRun(150,220000,0,0,i,"2").compute_hydrostatic_thickness()
    
    ht_tight = hd_tight[7]
    ht_wide = hd_wide[7]
    ht_ultra_wide = hd_ultra_wide[7]
    
    hd_list_tight.append(ht_tight)
    hd_list_wide.append(ht_wide)
    hd_list_ultra_wide.append(ht_ultra_wide)
    
    
    
plt.plot(hd_list_tight)
plt.plot(hd_list_wide)
plt.plot(hd_list_ultra_wide)