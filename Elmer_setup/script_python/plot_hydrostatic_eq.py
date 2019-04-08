#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 15:29:02 2019

@author: jloos
"""

import numpy as np
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
import turbulucid
from modelrun import ModelRun
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
import numpy as np
from tempfile import TemporaryFile
import numpy as np
import matplotlib
from pylab import *
import numpy as np
from pandas import DataFrame, Series
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
from scipy.interpolate import griddata,interp2d,interp1d
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay,ConvexHull
from matplotlib.tri import Triangulation,TriInterpolator, LinearTriInterpolator,CubicTriInterpolator, TriFinder
import pylab
from numpy.linalg import inv
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import LinearNDInterpolator,Rbf, UnivariateSpline,CloughTocher2DInterpolator, griddata, CubicSpline
import concave_hull
from collections import OrderedDict
from matplotlib.colors import LogNorm
from pylab import rcParams





class Plot_hydrostatic_deviation():
    
    """
    Class Plot:
    Here description
    """
    
# Check if volume is mounted, if not mount volume
    
    
    
    def __init__(self, 
                 **kwargs):
        
        
        """
        ModelRun 
        
        Parameters
        ----------
        bump_amplitude : int
            height of the bump
        bump_distribution_x : int
            height of the bump
        bump_distribution_y : int
            height of the bump    
        prop : string
            extra specification e.g. double (two gauss functions)
         timestep : int
            number of timestep
        """