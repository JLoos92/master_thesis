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


#import pykrige.kriging_tools as kt
###################### 0. Load and arrange data ##################################


#MeshResolution
meshresolution = 1000
#Choose model domain smaller than data domain with this value.
dist2boundary = 5000
#Domain Width in m+5km to avoide interpolation issues at boundary; spacing;
LY = 50 * 1000
dy = 1000; print (dy)

#Domain Length in m; Center at GL
LX = 250 * 1000
dx = 1000; print (dx)

#Grounding line x in Mismip2D
GLx = 1054 * 1000
#Vertical Layers for Velocity Interpolation
vertlayers = 21


#Load Input from Mismip2D
Mismip2D = genfromtxt('mismip_1a4_09.nh.csv', delimiter=',')

indnonzeroZs = np.where(Mismip2D[:,6] == 0)[0]
indnonzeroZb = np.where(Mismip2D[:,7] == 0)[0]
indnonzeroBed = np.where(Mismip2D[:,8] == 0)[0]


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
xv=np.arange(-LX/2,LX/2,dx)+GLx;
yv=np.arange(-LY/2,LY/2,dy);

#Bedi = interp1d(x(indnonzeroBed),Bed,xv,'PCHIP');
#Zbi = interp1d(x(indnonzeroZb),Zb,xv,'PCHIP');
Bed_i = interp1d(x[indnonzeroBed],Bed);
Zs_i = interp1d(x[indnonzeroZs],Zs);
Zb_i = interp1d(x[indnonzeroZb],Zb);

Zb_i_new = Zb_i(xv)
Zs_i_new = Zs_i(xv)
Bed_i_new = Bed_i(xv)
#numpy.interp(x, xp, fp, left=None, right=None)


#Get interpolated velocities at back boundary
#x2=x;
#[val_valindback] = min(abs(x2-min(xv)));
#backboundary1 = np.where(x2==x(valindback));
#x2(backboundary1)=NaN;
#[val2_valindback2] = min(abs(x2-min(xv)));
#backboundary2 = find(x2==x(valindback2));



#indnonzeroZs = find(Mismip2D(:,7)~=0);
#indnonzeroZb = find(Mismip2D(:,9)~=0);
#indnonzeroBed = find(Mismip2D(:,9)~=0);





