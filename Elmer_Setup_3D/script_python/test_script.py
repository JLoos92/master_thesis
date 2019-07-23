# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 15:55:13 2018
"""

widths = np.arange(10000,100000,10000)

org = []
for i in widths:
    
    
    original_width = int(round(np.sqrt(i)*2)*2)
    
    org.append(original_width)