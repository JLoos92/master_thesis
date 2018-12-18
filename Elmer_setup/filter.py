#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 10:15:05 2018

@author: Schmulius
"""

import numpy as np
import vtk
import mayavi
from mayavi import mlab


import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import sys

import matplotlib.pyplot as plt



# Try to read vtk files
gridreader = vtk.vtkXMLPUnstructuredGridReader()
gridreader.SetFileName("Forward5KM10L0151.pvtu")
gridreader.Update()
blabla = gridreader.GetOutputUpdateExtent ()
data = gridreader.GetOutput()
#gridreader.GetNumberOfPoints()
#gridreader.GetNumberOfCells()
#gridreader.GetNumberOfPointArrays()



#points = data.GetPoints()
#npts = points.GetNumberOfPoints()
#x = vtk_to_numpy(points.GetData())

#triangles=  vtk_to_numpy(data.GetCells().GetData())
#triangles[4:8]
#ntri = triangles.size//4  # number of cells
#tri = np.take(triangles,[n for n in range(triangles.size) if n%4 != 0]).reshape(ntri,3)



# saved by elmer



#plt.rcParams['figure.figsize'] = (10.0, 6.0)
#n_arrays = reader.GetNumberOfPointArrays()
#for i in range(n_arrays):
 #   print(reader.GetPointArrayName(i))
 
    
#u = vtk_to_numpy(data.GetPointData().GetArray(2))



#plt.figure(figsize=(8, 8))
#plt.triplot(x[:,0], x[:,1], tri)
#plt.gca().set_aspect('equal')


