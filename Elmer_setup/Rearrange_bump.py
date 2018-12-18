#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 15:39:36 2018

@author: Schmulius
"""

###################### Packages ##################################

import matplotlib
#matplotlib.use('nbagg')
from pylab import *
import numpy as np
from pandas import DataFrame, Series
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
import plotly as py
py.tools.set_credentials_file(username='Julius92Loos', api_key='ZdBblhNTi4fwVZUtaGwy')
import plotly.graph_objs as go

# more
from scipy.stats import spearmanr
from scipy.stats import norm 
from scipy.stats import gamma
from scipy.stats import beta
from scipy.spatial import distance_matrix
from scipy import stats as sst
from scipy.stats import multivariate_normal
from scipy.spatial.distance import cdist
from IPython.display import Image
from astropy.io import ascii
from IPython.display import display
from IPython.display import SVG
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import axes3d, Axes3D #<-- Note the capitalization! 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter



###################### 0. Load and arrange data ##################################


ZB = genfromtxt('Mismip3DSetUpSteadyState_Remesh/DEM/ZB.xyz')
BED = genfromtxt('Mismip3DSetUpSteadyState_Remesh/DEM/BED.xyz')
ZS =  genfromtxt('Mismip3DSetUpSteadyState_Remesh/DEM/ZS.xyz')
 
# Rearrange ZB   
y= np.unique(ZB[:,0])
x= np.unique(ZB[:,1])

xx,yy = np.meshgrid(x, y)

ZB_B = np.asmatrix(ZB[:,2])
ZB_B = ZB_B.reshape(201,51)
    

# Rearrange BED

BED_B = np.asmatrix(BED[:,2])
BED_B = BED_B.reshape(201,51)

# Rearrange ZS

ZS_B = np.asmatrix(ZS[:,2])
ZS_B = ZS_B.reshape(201,51)


 
# Setup for Bumpfunction
GLx = 1055 * 1000
# Add Gauss Bump to BED topography

# Bump function with sigma and amplitude

M = np.zeros((201,51))
N = 1; # number of bumps
maxAmplitude = 60
sigmax = 200
sigmay = 500
theta = 2*np.pi     #rotation of bump

dl = 0

a = ((np.cos(theta))**2)/(2*(sigmax**2)) + ((np.sin(theta))**2)/(2*sigmay**2)
b = (-np.sin(2*theta))/(4*(sigmax**2)) + (np.sin(2*theta))/(4*sigmay**2)
c = ((np.sin(theta))**2)/(2*(sigmax**2)) + ((np.cos(theta))**2)/(2*(sigmay**2))

  
#location of bumps or xc, yc for grounding line position
xc = 0
dl = 0
yc = GLx - dl
    




#Gauss function with theta and coeffecient for rotation
exponent_el = (((a*(xx-xc)**2)+ (2*b*(xx-xc)*(yy-yc)) + (c*(yy-yc)**2)))
amplitude = maxAmplitude
    
 # add Gauss to the matrix M
 
Zero = np.zeros((201,51)) 
 

M = M + (amplitude*np.exp(-exponent_el))


# Add Gauss to topography (plot matrice for regular and rot)
BED_B_PLOT = BED_B + M
ZB_B_PLOT = ZB_B + M




# Arrange Data for ELMER/ICE


# Cut ZB for ELMER

M_ZB = ZB_B_PLOT[0:101,:]
M_Zero = ZB_B[101:201,:]

ZB_B_PLOT = np.append(M_ZB,M_Zero,axis=0)


#from np.matrix to np.array
BED_B_safe = np.squeeze(np.asarray(BED_B_PLOT))
ZB_B_safe = np.squeeze(np.asarray(ZB_B_PLOT))

BED_B_safe = BED_B_safe.flatten()
ZB_B_safe = ZB_B_safe.flatten()

ZB_B_safe = (np.matrix([BED[:,0],BED[:,1],ZB_B_safe])).T
BED_B_safe = (np.matrix([BED[:,0],BED[:,1],BED_B_safe])).T


# Write new DEM file for BED surface
np.savetxt("BED_bump" + str(maxAmplitude)+ "_" + str(sigmax) + str(sigmay)+ "_" + str(dl) + ".xyz",BED_B_safe)
np.savetxt("ZB_bump" + str(maxAmplitude)+ "_" + str(sigmax) + str(sigmay)+ "_" + str(dl) + ".xyz",ZB_B_safe)



   
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(xx,yy,BED_B_PLOT,cmap=cm.coolwarm)
ax.plot_surface(xx,yy,ZS_B,cmap=cm.coolwarm)
ax.plot_surface(xx,yy,ZB_B_PLOT,cmap=cm.coolwarm)#ax.plot_surface(xx,yy,ZS_B,cmap=cm.coolwarm)
ax.view_init(azim=-151,elev=23)
   
