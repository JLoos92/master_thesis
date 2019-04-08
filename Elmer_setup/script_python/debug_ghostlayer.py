#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 16:11:30 2019

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
from collections import OrderedDict
from matplotlib.colors import LogNorm
from pylab import rcParams
import vtk
from vtk.util.numpy_support import vtk_to_numpy

Dict = {
'ghostlayer' : ModelRun(250,300,300,'double',50).compute_hydrostatic_thickness()
        }

O_Dict = OrderedDict(Dict) 




# Define properties for subplots
    
font_title = {'color':'black',
       'size':'40',
       'weight':'bold'
            }

font_axes = {'color':'black',
       'size':'14',
       'weight':'italic'
           }

font_annotation = {'color':'black',
                   'size': '10'
                   }

# Make subplots and iterate over dictionary for hydrostatic imbalances

run = 'ghostlayer'
x = O_Dict[run][0]
y = O_Dict[run][1]
ht =O_Dict[run][2]
   
plt.tripcolor(x,y,ht,shading='gouraud',vmin=-15,vmax=15,cmap = 'RdBu')   
plt.show()


cutter = ModelRun(250,300,300,'double',50).cutter()

velocity = cutter.GetOutput().GetPointData().GetArray('dsdt')
velocity = vtk_to_numpy(velocity)
points =  cutter.GetOutput().GetPoints().GetData()
points = vtk_to_numpy(points)



x1 = points[:,0]
y1 = points[:,1]
vel = velocity
   
#plt.tripcolor(x1,y1,vel,shading='gouraud',cmap = 'RdBu') 
#plt.show()  

f_name = '/Volumes/esd01/docs/jloos/data_small/runs_elmerice_fixed/Mesh250_200200_double/Mesh/ForwardGLFixed0015.pvtu'
xmlReader = vtk.vtkXMLPUnstructuredGridReader()
xmlReader.SetFileName(f_name)
xmlReader.Update()


def selectMaxz(x, y, z):
            # Get grouped indices
            sidx = (y + x*(y.max() - y.min() + 1)).argsort()
           
        
            # sort x, y, z
            x_sorted = x[sidx]
            y_sorted = y[sidx]
            z_sorted = z[sidx]
        
            # Get equality mask between each sorted X and Y elem against previous 
            seq_eq_mask = (x_sorted[1:] == x_sorted[:-1]) & (y_sorted[1:] == y_sorted[:-1])
            cut_idx = np.flatnonzero(np.concatenate(( [True], ~seq_eq_mask)))
        
            # Use those cut_idx to get intervalled maximum values
            maxZ = np.minimum.reduceat(z_sorted, cut_idx)
        
            # Make tuples of the groupings of x,y and the corresponding min Z values
            return (zip(x_sorted[cut_idx], y_sorted[cut_idx]), maxZ.tolist())
        
        
max_points = selectMaxz(points[:,0],points[:,1],points[:,2])       
z_max = max_points[1]
z = points[:,2]
zmax_ind = np.isin(z,z_max)
zmax_ind_num = np.where(zmax_ind)
vel_new = vel[zmax_ind]
x_new = x1[zmax_ind]
y_new = y1[zmax_ind]


plt.tripcolor(x_new,y_new,vel_new,shading='gouraud',cmap = 'RdBu') 
plt.show() 


x =  points[:,0]
y =  points[:,1]
z =  points[:,2]

p_w = 1000.0 #kg m−3 ), ice (ρi =918kgm−3), and air (ρa =2kgm−3):
p_i = 900.0

zs = vtk_to_numpy(cutter.GetOutput().GetPointData().GetArray('zs'))
zb = vtk_to_numpy(cutter.GetOutput().GetPointData().GetArray('zb'))

def selectMinz(x, y, z):
            # Get grouped lex-sort indices
            sidx = (y + x*(y.max() - y.min() + 1)).argsort()
            
        
            # Sort x, y, z
            x_sorted = x[sidx]
            y_sorted = y[sidx]
            z_sorted = z[sidx]
        
            # Get equality mask between each sorted X and Y elem against previous 
            seq_eq_mask = (x_sorted[1:] == x_sorted[:-1]) & (y_sorted[1:] == y_sorted[:-1])
            cut_idx = np.flatnonzero(np.concatenate(( [True], ~seq_eq_mask)))
        
            # Use those cut_idx to get intervalled minimum values
            minZ = np.minimum.reduceat(z_sorted, cut_idx)
        
            # Make tuples of the groupings of x,y and the corresponding min Z values
            return ((x_sorted[cut_idx], y_sorted[cut_idx]), minZ.tolist())
        
         # Function to find min-values of z_s
         
def selectMaxz(x, y, z):
            # Get grouped indices
            sidx = (y + x*(y.max() - y.min() + 1)).argsort()
           
        
            # sort x, y, z
            x_sorted = x[sidx]
            y_sorted = y[sidx]
            z_sorted = z[sidx]
        
            # Get equality mask between each sorted X and Y elem against previous 
            seq_eq_mask = (x_sorted[1:] == x_sorted[:-1]) & (y_sorted[1:] == y_sorted[:-1])
            cut_idx = np.flatnonzero(np.concatenate(( [True], ~seq_eq_mask)))
        
            # Use those cut_idx to get intervalled maximum values
            maxZ = np.maximum.reduceat(z_sorted, cut_idx)
        
            # Make tuples of the groupings of x,y and the corresponding min Z values
            return (zip(x_sorted[cut_idx], y_sorted[cut_idx]), maxZ.tolist())
        

zmin = selectMinz(x,y,z)
min_pairs = selectMinz(x,y,z)
min_pairs = min_pairs[0]
#min_pairs = np.asanyarray(min_pairs)


x_corr = min_pairs[0] 
y_corr = min_pairs[1]
zmax = selectMaxz(x,y,z)

zmax = zmax[1]
zmax = np.asarray(zmax)

zmin = zmin[1]
zmin = np.asarray(zmin)


zmin_ind = np.isin(z,zmin)
zmax_ind = np.isin(z,zmax)

zmin_ind = np.where(zmin_ind)
zmax_ind = np.where(zmax_ind)

zb_new = zb[zmin_ind]
zs_new = zs[zmax_ind]

zb_new = zmin
zs_new = zmax
#x_new = x[zmax_ind]
#y_new = y[zmax_ind]

thick_model = -zs_new+zb_new       
thick_calc = np.divide((p_w*zs_new),(p_w-p_i))                
h_thickness = thick_model + thick_calc 



plt.tripcolor(x_corr,y_corr,h_thickness,shading='gouraud',vmin=-30,vmax=30,cmap = 'RdBu') 
plt.show() 

plt.tripcolor(x_new,y_new,zb_new,shading='gouraud',cmap = 'RdBu')
plt.title('old function') 
plt.show() 


plt.tripcolor(x_corr,y_corr,zmin,shading='gouraud',cmap = 'RdBu') 
plt.title('ZB')
plt.show()

plt.tripcolor(x_new,y_new,zs_new,shading='gouraud',cmap = 'RdBu') 
plt.title('ZS')
plt.show()

plt.tripcolor(x_corr,y_corr,zmax,shading='gouraud',cmap = 'RdBu') 
plt.title('ZS_real')
plt.show()
