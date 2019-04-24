#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 14:18:43 2019

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
from scipy.spatial import Delaunay,ConvexHull
from matplotlib.tri import Triangulation,TriInterpolator, LinearTriInterpolator,CubicTriInterpolator, TriFinder
import pylab
from numpy.linalg import inv
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import LinearNDInterpolator,Rbf, UnivariateSpline,CloughTocher2DInterpolator
import concave_hull

matplotlib.get_backend()
matplotlib.rcParams['figure.figsize'] = (25,15)


#cutter = ModelRun(50,500,500,0,0).cutter()
convexhull = ModelRun(50,500,500,0,50).convexhull()
cutter = ModelRun(50,500,500,0,50).cutter()
dictvar = ModelRun(250,500,500,0,100).dict_var

points = vtk_to_numpy(convexhull.GetPoints().GetData())
#_cutter = vtk_to_numpy(cutter.GetPoints().GetData())
zs = vtk_to_numpy(convexhull.GetPointData().GetArray('zs'))
zb = vtk_to_numpy(convexhull.GetPointData().GetArray('zb'))
depth = vtk_to_numpy(convexhull.GetPointData().GetArray('depth'))
#stress = vtk_to_numpy(convexhull.GetPointData().GetArray('stress vector'))
zs = ma.masked_where(zs<=0,zs)
x = points[:,0]
y = points[:,1]
z = points[:,2]




fig = plt.figure()
ax = Axes3D(fig)

ax.scatter(x,y,zs)
plt.show()


tri = Triangulation(x,y)


tri = Triangulation(x,y)
points = np.delete(points, 2, 1)
delaunay = Delaunay(points)

#------------------------------------------------------------------------------
#   Line Profiles
#------------------------------------------------------------------------------

#Cutter1
line=vtk.vtkPlane()
line.SetOrigin(1065000,0,0)
line.SetNormal(1,0,0)


cutter1 = vtk.vtkCutter()   
cutter1.SetCutFunction(line)
cutter1.SetInputConnection(cutter.GetOutputPort())
cutter1.Update()
line_points = vtk_to_numpy(cutter1.GetOutput().GetPoints().GetData())

#Cutter2
line2=vtk.vtkPlane()
line2.SetOrigin(1070000,0,0)
line2.SetNormal(1,0,0)


cutter2 = vtk.vtkCutter()   
cutter2.SetCutFunction(line2)
cutter2.SetInputConnection(cutter.GetOutputPort())
cutter2.Update()
line_points2 = vtk_to_numpy(cutter2.GetOutput().GetPoints().GetData())


#Cutter3
line3=vtk.vtkPlane()
line3.SetOrigin(1075000,0,0)
line3.SetNormal(1,0,0)


cutter3 = vtk.vtkCutter()   
cutter3.SetCutFunction(line3)
cutter3.SetInputConnection(cutter.GetOutputPort())
cutter3.Update()
line_points3 = vtk_to_numpy(cutter3.GetOutput().GetPoints().GetData())
stress = vtk_to_numpy(cutter3.GetOutput().GetPointData().GetArray('stress vector'))
stress_z = stress[:,2]
#------------------------------------------------------------------------------
#   Convex Hull for crosssection
#------------------------------------------------------------------------------
hull1_points = np.delete(line_points, 0, 1)
hull1 = concave_hull.compute(hull1_points,2)

hull2_points = np.delete(line_points2, 0, 1)
hull2 = concave_hull.compute(hull2_points,3)

hull3_points = np.delete(line_points3, 0, 1)
hull3 = concave_hull.compute(hull3_points,4)



font ={'color':'red',
       'size':'20'
       }



#plt.tripcolor(tri,ctd_values,shading='gouraud')
plt.plot(line_points[:,0],line_points[:,1],'-r',label = 'A')
plt.text(1065000,-5000,'A',fontdict=font)
plt.text(1065000,5000,"A'",fontdict=font)
plt.plot(line_points2[:,0],line_points2[:,1],'-r',label = 'B')
plt.text(1070000,-5000,'B',fontdict=font)
plt.text(1070000,5000,"B'",fontdict=font)
plt.plot(line_points3[:,0],line_points3[:,1],'-r',label = 'C')
plt.text(1075000,-5000,'C',fontdict=font)
plt.text(1075000,5000,"C'",fontdict=font)






plt.triplot(x,y,delaunay.simplices,alpha=0.2)
plt.tripcolor(tri,zs,shading='gouraud')
plt.title('Delaunay triangulation of refined mesh')
plt.clim(30,35)
plt.colorbar()
plt.gca().set_aspect('equal')
plt.show()



p_w = 1027 #kg m−3 ), ice (ρi =918kgm−3), and air (ρa =2kgm−3):
p_i = 900
p_a = 2
H_a = 0


zcalc = np.asarray([(i>0) * i for i in z])
hh = zs - zcalc 
hh = np.asarray([(i>0) * i for i in hh])
hh = np.divide((p_w*zs),(p_w-p_i))



# Profiles of channel imprint

#plt.scatter(line_points[:,1],line_points[:,2],label = 'A')
#plt.scatter(line_points2[:,1],line_points2[:,2],label = 'B')
#plt.scatter(line_points3[:,1],line_points3[:,2],label = 'C')
#plt.plot(hull1_points[hull1.vertices,0],hull1_points[hull1.vertices,1],'r--',lw=2)
#plt.legend(loc='upper left')


#plt.plot(line_points[:,1],line_points[:,2],'o')

#for simplex in hull1.simplices:
        
 #   plt.plot(hull1_points[simplex,0],hull1_points[simplex,1],'k-')
    
#plt.plot(hull1_points[hull1.vertices,0],hull1_points[hull1.vertices,1],'r--',lw=2)




plt.plot(hull1[:,0],hull1[:,1],'r-')
plt.plot(hull2[:,0],hull2[:,1],'b-')
plt.plot(hull3[:,0],hull3[:,1],'g-')
plt.plot(line_points3[:,1],line_points3[:,2],stress_z)
plt.title('Channel imprints of crosssection')
plt.grid()
plt.show()
