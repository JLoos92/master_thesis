#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:29:50 2019

@author: jloos
"""

# rms 
#def rmsvalue(arr,n):
#    square = 0
#    mean = 0.0
#    root= 0.0
#    
#    for i in test(0,n):
#        square = (arr[i]**2)
#        print(square)
#        mean = (square / (float)(n))
#        print('mean = ', mean)
#        root = math.sqrt(mean)
#        
#        return root

test = np.linspace(1,20,10)
n = len(test)
#value = rmsvalue(test,n)


rms = np.sqrt(np.mean(test**2))