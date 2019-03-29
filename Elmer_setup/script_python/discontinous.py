#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 16:10:11 2019

@author: jloos
"""

from modelrun import *
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
from scipy.interpolate import BSpline,splrep,interp1d, interp2d
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from matplotlib.tri import Triangulation,TriInterpolator, LinearTriInterpolator,CubicTriInterpolator, TriFinder
import pylab
from numpy.linalg import inv
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import LinearNDInterpolator,Rbf
from scipy.ndimage.filters import gaussian_filter

matplotlib.rcParams['figure.figsize'] = (25,15)



c_hull = ModelRun(150,500,500,0,150).compute_concavehull(1076000,3)
cal = ModelRun(150,500,500,0,150).cut_and_slice(1076000,'stress vector')

stress = cal[0]
linepoints = cal[1]

x_smooth = np.linspace(c_hull[:,0].min(),c_hull[:,0].max(), 100)
x = c_hull[:,0]
y = c_hull[:,1]

xx = linepoints[:,1]
yy = linepoints[:,2]

plt.plot(x,y,)
plt.scatter(linepoints[:,1],linepoints[:,2],s = stress[:,0]/10000)
plt.grid()

X,Y = np.meshgrid(xx,yy)

f = interp2d(xx,yy,stress[:,0])


z = gaussian_filter(stress[:,0],sigma=5,mode='wrap')

plt.pcolormesh(xx,yy,z)






