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
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from matplotlib.tri import Triangulation,TriInterpolator, LinearTriInterpolator,CubicTriInterpolator, TriFinder
import pylab
from numpy.linalg import inv
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import LinearNDInterpolator,Rbf

matplotlib.rcParams['figure.figsize'] = (25,15)











cutter = ModelRun(250,500,500,0,100).cutter()
cutter_1 = ModelRun(250,500,500,0,1).cutter()
run = ModelRun(250,500,500,0,100).xmlReader




#xmlReader = vtk.vtkXMLPUnstructuredGridReader()
#xmlReader.SetFileName(f_name)
#xmlReader.Update()


points_cut = vtk_to_numpy(cutter.GetOutput().GetPoints().GetData())
points_cut_1 = vtk_to_numpy(cutter_1.GetOutput().GetPoints().GetData())
velo_cut= vtk_to_numpy(cutter.GetOutput().GetPointData().GetArray(17))
velo_cut= velo_cut[:,0]

zs = vtk_to_numpy(cutter.GetOutput().GetPointData().GetArray('zs'))
dsdt = vtk_to_numpy(cutter.GetOutput().GetPointData().GetArray('dsdt'))




# Extract coordinates
x = points_cut[:,0]
y = points_cut[:,1]
z = points_cut[:,2]

# Extract coordinates
x_1 = points_cut_1[:,0]
y_1 = points_cut_1[:,1]
z_1 = points_cut_1[:,2]


tri = Triangulation(x,y)
tri_1 = Triangulation(x_1,y_1)
points = np.delete(points_cut, 2, 1)
points_1 = np.delete(points_cut_1, 2, 1)  
delaunay = Delaunay(points)
delaunay_1 = Delaunay(points_1)
#
plt.triplot(points_cut[:,0],points_cut[:,1],delaunay.simplices,alpha=0.2)
plt.tripcolor(tri,dsdt,shading='gouraud')
plt.title('Delaunay triangulation of refined mesh')
#plt.clim(20,40)
plt.colorbar()
plt.gca().set_aspect('equal')
plt.show()



p_w = 1027 #kg m−3 ), ice (ρi =918kgm−3), and air (ρa =2kgm−3):
p_i = 900
p_a = 2
H_a = 0

# Parameters for calculating hydrostatic thickness

h_100_calculated = ((p_w*zs)/(p_w-p_i))-(((p_a-p_i)/(p_i-p_w))* H_a)
h_001_calculated  = ((p_w*z_1)/(p_w-p_i))-(((p_a-p_i)/(p_i-p_w))*H_a)


#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#Cutter1
line=vtk.vtkPlane()
line.SetOrigin(1065000,0,0)
line.SetNormal(1,0,0)


cutter = vtk.vtkCutter()   
cutter.SetCutFunction(line)
cutter.SetInputConnection(run.GetOutputPort())
cutter.Update()
line_points = vtk_to_numpy(cutter.GetOutput().GetPoints().GetData())

#Cutter2
line2=vtk.vtkPlane()
line2.SetOrigin(1070000,0,0)
line2.SetNormal(1,0,0)


cutter2 = vtk.vtkCutter()   
cutter2.SetCutFunction(line2)
cutter2.SetInputConnection(run.GetOutputPort())
cutter2.Update()
line_points2 = vtk_to_numpy(cutter2.GetOutput().GetPoints().GetData())


#Cutter3
line3=vtk.vtkPlane()
line3.SetOrigin(1075000,0,0)
line3.SetNormal(1,0,0)


cutter3 = vtk.vtkCutter()   
cutter3.SetCutFunction(line3)
cutter3.SetInputConnection(run.GetOutputPort())
cutter3.Update()
line_points3 = vtk_to_numpy(cutter3.GetOutput().GetPoints().GetData())



fig = plt.figure()
plt.triplot(points_cut[:,0],points_cut[:,1],delaunay.simplices,alpha=0.2)
plt.tripcolor(tri,dsdt,shading='gouraud')
plt.plot(line_points[:,0],line_points[:,1],'-r')
plt.text(1065000,-5000,'A')
plt.text(1065000,5000,'A'')
plt.plot(line_points2[:,0],line_points2[:,1],'-r')
plt.plot(line_points3[:,0],line_points3[:,1],'-r')
plt.text(1070000,-5000,'B')
plt.text(107000,5000,'B'')
plt.text(1075000,-5000,'C')
plt.text(107500,5000,'C'')
plt.legend()

plt.title('Delaunay triangulation of refined mesh')
#plt.clim(20,40)
plt.colorbar()
plt.gca().set_aspect('equal')
plt.show()




#fig = plt.figure()
#plt.triplot(points_cut[:,0],points_cut[:,1],delaunay.simplices,alpha=0.1)
#plt.tripcolor(tri,h_100_calculated,shading='gouraud')
#plt.plot(line_points[:,0],line_points[:,1],'-r')
#plt.text(1065000,-5000,'A')
#plt.plot(line_points2[:,0],line_points2[:,1],'-r')
#plt.plot(line_points3[:,0],line_points3[:,1],'-r')
#plt.legend()
#
#plt.title('Delaunay triangulation of refined mesh')
##plt.clim(20,40)
#plt.colorbar()
#plt.gca().set_aspect('equal')
#plt.show()

#plt.triplot(points_cut_1[:,0],points_cut_1[:,1],delaunay_1.simplices,alpha=0.1)
#plt.tripcolor(tri_1,h_001_calculated,shading='gouraud')
#plt.title('Delaunay triangulation of refined mesh')
##plt.clim(20,40)
#plt.colorbar()
#plt.gca().set_aspect('equal')
#plt.show()










