
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:36:40 2019

@author: Schmulius
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
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from matplotlib.tri import Triangulation,TriInterpolator, LinearTriInterpolator,CubicTriInterpolator, TriFinder
import pylab
from numpy.linalg import inv
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import LinearNDInterpolator,Rbf



cutter = ModelRun(250,250,250,0,50).cutter()
dictvar = ModelRun(250,250,250,0,50).dict_var

points = vtk_to_numpy(cutter.GetOutput().GetPoints().GetData())
x = points[:,0]
y = points[:,1]
z = points[:,2]

def selectMinz_vectorized(x, y, z):
    # Get grouped lex-sort indices
    sidx = (y + x*(y.max() - y.min() + 1)).argsort()
    # or sidx = np.lexsort([x, y])

    # Lex-sort x, y, z
    x_sorted = x[sidx]
    y_sorted = y[sidx]
    z_sorted = z[sidx]

    # Get equality mask between each sorted X and Y elem against previous ones.
    # The non-zero indices of its inverted mask gives us the indices where the 
    # new groupings start. We are calling those as cut_idx.
    seq_eq_mask = (x_sorted[1:] == x_sorted[:-1]) & (y_sorted[1:] == y_sorted[:-1])
    cut_idx = np.flatnonzero(np.concatenate(( [True], ~seq_eq_mask)))

    # Use those cut_idx to get intervalled minimum values
    minZ = np.minimum.reduceat(z_sorted, cut_idx)

    # Make tuples of the groupings of x,y and the corresponding min Z values
    return (zip(x_sorted[cut_idx], y_sorted[cut_idx]), minZ.tolist())


def selectMaxz_vectorized(x, y, z):
    # Get grouped lex-sort indices
    sidx = (y + x*(y.max() - y.min() + 1)).argsort()
    # or sidx = np.lexsort([x, y])

    # Lex-sort x, y, z
    x_sorted = x[sidx]
    y_sorted = y[sidx]
    z_sorted = z[sidx]

    # Get equality mask between each sorted X and Y elem against previous ones.
    # The non-zero indices of its inverted mask gives us the indices where the 
    # new groupings start. We are calling those as cut_idx.
    seq_eq_mask = (x_sorted[1:] == x_sorted[:-1]) & (y_sorted[1:] == y_sorted[:-1])
    cut_idx = np.flatnonzero(np.concatenate(( [True], ~seq_eq_mask)))

    # Use those cut_idx to get intervalled maximum values
    maxZ = np.maximum.reduceat(z_sorted, cut_idx)

    # Make tuples of the groupings of x,y and the corresponding min Z values
    return (zip(x_sorted[cut_idx], y_sorted[cut_idx]), maxZ.tolist())



zmin = selectMinz_vectorized(x,y,z)
zmax = selectMaxz_vectorized(x,y,z)

zmax = zmax[1]
zmax = np.asarray(zmax)

zmin = zmin[1]
zmin = np.asarray(zmin)


zs = vtk_to_numpy(cutter.GetOutput().GetPointData().GetArray('zs'))
zb = vtk_to_numpy(cutter.GetOutput().GetPointData().GetArray('zb'))


zmin_ind = np.isin(z,zmin)
zmax_ind = np.isin(z,zmax)

zmin_ind = np.where(zmin_ind)
zmax_ind = np.where(zmax_ind)

zb_new = zb[zmin_ind]
zs_new = zs[zmax_ind]

zb_s = zb[zb!=0]
zs_s = zs[zs!=0]