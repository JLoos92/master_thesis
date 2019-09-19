#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 11:25:08 2018

@author: Julius Loos
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
print ("##################################### \nThis script will compile BED.xyz and ZB.xyz for ELMER/ICE with a given gaussian bump amplitude and a specified sigma.\
 Choose between 140 - 900 for amplitude and 1000 - 10000 for distribution \n##################################### ")

#ax = plt.axes(projection='3d')
#MeshResolution
meshresolution = 1000;
#Choose model domain smaller than data domain with this value.
dist2boundary = 5000;
#Domain Width in m+5km to avoide interpolation issues at boundary; spacing;
LY = 51 * 1000;
dy = 1000

#Domain Length in m; Center at GL
LX = 251 * 1000;
dx = 1000

#Grounding line x in Mismip2D
GLx = 1054 * 1000
#Vertical Layers for Velocity Interpolation
vertlayers = 21


#Load Input from Mismip2D
Mismip2D = genfromtxt('mismip_1a4_09.nh.csv', delimiter=',')

indnonzeroZs = np.where(Mismip2D[:,6] != 0)[0]
indnonzeroZb = np.where(Mismip2D[:,8] != 0)[0]
indnonzeroBed = np.where(Mismip2D[:,8] != 0)[0]


Zs = Mismip2D[indnonzeroZs,6]
V1 = Mismip2D[:,19];
V2 = Mismip2D[:,20];
Depth = Mismip2D[:,10];
Height = Mismip2D[:,11];
Zb = Mismip2D[indnonzeroZb,8];
Bed = Mismip2D[indnonzeroBed,5];
x = Mismip2D[:,25];
z = Mismip2D[:,26];



#Get interpolated  2D DEM centered at GL and clipped to area of interest
xv=np.arange(-LX/2,LX/2,dx)+GLx+500;
xv=np.array(xv)
xv=np.array([xv])
yv=np.arange(-LY/2,LY/2,dy)+500;
yv=np.array(yv)
yv=np.array([yv])

hallo = np.arange(-LY/2,LY/2,dy)


Bed_i = interp1d(x[indnonzeroBed],Bed);
Zs_i = interp1d(x[indnonzeroZs],Zs);
Zb_i = interp1d(x[indnonzeroZb],Zb);

Zb_i_new = Zb_i(xv)
Zs_i_new = Zs_i(xv)
Bed_i_new = Bed_i(xv)

Zb = np.array([Zb])

#Get interpolated velocities at back boundary
x2=x;
val = np.min(np.abs(x2 - np.min(xv)))
valindex = np.abs(x2 - np.min(xv));
valindex = np.argmin(valindex)
backboundary1 = np.where(x2==x[valindex]);
x2[backboundary1] = float('nan');


val2 = np.nanmin(np.abs(x2 - np.min(xv)))
valindex2 = np.abs(x2 - np.min(xv));
valindex2 = np.nanargmin(valindex2)
backboundary2 = np.where(x2==x[valindex2]);


meanV1back = (val*V1[backboundary1]+val2*V1[backboundary2])/(val+val2);
meanV2back = (val*V2[backboundary1]+val2*V2[backboundary2])/(val+val2);
meanHeightback = (val*Height[backboundary1]+val2*Height[backboundary2])/(val+val2);

# Fit Polynomial, to be used at back boundary in Elmer
coeffV1 = np.polyfit(meanHeightback,meanV1back,4);
coeffV2 = np.polyfit(meanHeightback,meanV2back,2);



#Extrusion
xx,yy = np.meshgrid(yv,xv)
xx# = xx.T
yy #= yy.T

ZB = np.zeros((251,yv.size)) 
BED = np.zeros((251,yv.size))  
ZS = np.zeros((251,yv.size))  
 
for k in range(yv.size):

    ZB[:,k] = Zb_i_new
    BED[:,k] = Bed_i_new
    ZS[:,k] = Zs_i_new

# Write DEM Files to be read into Elmer
  


ZBout = ZB.flatten()
BEDout = BED.flatten()
ZSout = ZS.flatten()
xxout = xx.flatten()
yyout = yy.flatten()

ZBout = (np.matrix([yyout,xxout,ZBout])).T
BEDout = (np.matrix([yyout,xxout,BEDout])).T
ZSout = (np.matrix([yyout,xxout,ZSout])).T

#from np.matrix to np.array
ZBout = np.squeeze(np.asarray(ZBout))
BEDout = np.squeeze(np.asarray(BEDout))
ZSout = np.squeeze(np.asarray(ZSout))

#Sorting
#ZBouts = ZBout[ZBout[:,0].argsort(),]
#BEDouts = BEDout[BEDout[:,0].argsort(),]
#ZSouts = ZSout[ZSout[:,0].argsort(),]


# Plots
#plt.plot(x,z,'r.', label='Original data')

#plt.ylabel('Height (m)')
#plt.xlabel('Distance (km)')
#plt.show()

###################### 1. Setup for bump ##################################

# Add Gauss Bump to BED topography

# Bump function with sigma and amplitude

M = np.zeros((251,51))
N = 1; # number of bumps
sigma = 500
maxAmplitude = 200
sigmax = 500
sigmay = 1000
theta = 2*np.pi     #rotation of bump

#dl = 400

a = ((np.cos(theta))**2)/(2*(sigmax**2)) + ((np.sin(theta))**2)/(2*sigmay**2)
b = (-np.sin(2*theta))/(4*(sigmax**2)) + (np.sin(2*theta))/(4*sigmay**2)
c = ((np.sin(theta))**2)/(2*(sigmax**2)) + ((np.cos(theta))**2)/(2*(sigmay**2))

  
#location of bumps or xc, yc for grounding line position
xc = 0
dl = 0
yc = GLx - dl
    



# Gauss function normal
exponent = (((xx-xc)**2)+((yy-yc)**2))/(2*(sigma**2))

#Gauss function with theta and coeffecient for rotation
exponent_el = (((a*(xx-xc)**2)+ (2*b*(xx-xc)*(yy-yc)) + (c*(yy-yc)**2)))
amplitude = maxAmplitude
    
 # add Gauss to the matrix M
 
Zero = np.zeros((251,51)) 
 
M = M + (amplitude*np.exp(-exponent))
M_rot = M + (amplitude*np.exp(-exponent_el))


# Add Gauss to topography (plot matrice for regular and rot)
BED_B_PLOT = BED + M
ZB_B_PLOT = ZB + M

BED_B_PLOT_rot = BED + M_rot

# Cut ZB for ELMER

M_ZB = M_rot[0:125,:]
M_Zero = Zero[125:251,:]

ZB_B_PLOT_rot = ZB  + np.append(M_ZB,M_Zero,axis=0)






# Arrange Data for ELMER/ICE
BED_B = np.matrix.flatten(BED_B_PLOT)
BED_B = (np.matrix([yyout,xxout,BED_B])).T

BED_B_rot = np.matrix.flatten(BED_B_PLOT_rot)
BED_B_rot = (np.matrix([yyout,xxout,BED_B_rot])).T

ZB_B = np.matrix.flatten(ZB_B_PLOT)
ZB_B = (np.matrix([yyout,xxout,ZB_B])).T

ZB_B_rot = np.matrix.flatten(ZB_B_PLOT_rot)
ZB_B_rot = (np.matrix([yyout,xxout,ZB_B_rot])).T

#from np.matrix to np.array
BED_B = np.squeeze(np.asarray(BED_B))
ZB_B = np.squeeze(np.asarray(ZB_B))

BED_B_rot = np.squeeze(np.asarray(BED_B_rot))
ZB_B_rot = np.squeeze(np.asarray(ZB_B_rot))






# Write new DEM file for BED surface
#np.savetxt("BED_bump" + str(maxAmplitude) + str(sigmax) + str(sigmay)+ "_" + str(dl) + ".xyz",BED_B_rot)
#np.savetxt("ZB_bump" + str(maxAmplitude) + str(sigmax) + str(sigmay)+ "_" + str(dl) + ".xyz",ZB_B_rot)




###################### 2. Plots surface + bump ##################################


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(xx,yy,BED_B_PLOT_rot,cmap=cm.coolwarm)
ax.plot_surface(xx,yy,ZS,cmap=cm.coolwarm)
ax.plot_surface(xx,yy,ZB_B_PLOT_rot,cmap=cm.coolwarm)
ax.view_init(azim=172,elev=17)
