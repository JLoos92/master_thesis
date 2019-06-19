#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:33:32 2019

@author: jloos
"""

from __main__ import ModelRun
import pandas as pd
from scipy.interpolate import griddata

def printkwargs(**kwargs):
    
    
    kwargs = {
            'hi':2}



    print(kwargs)