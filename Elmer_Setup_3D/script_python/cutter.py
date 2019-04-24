#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 21:08:18 2019

@author: Schmulius
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import vtk
from vtk.util.misc import vtkGetDataRoot

import numpy as np
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.numpy_support import numpy_to_vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.interpolate import interp1d
from collections import OrderedDict

p1 = (50,500)
p2 = (50,0)

# Define points p1 and p2 as min,max for line-plot
p1 = np.append(np.array(p1), xmlReader2.vtkData.GetPoints()[0, 2])
p2 = np.append(np.array(p2), 0)

   # Compute the plane-normal as a cross-product
unit = (p2 - p1)/np.linalg.norm(p2 - p1)
geomNormal = np.array([0, 0, 1])
planeNormal = np.cross(unit, geomNormal)

    # Define the cutting plane
plane = vtk.vtkPlane()
plane.SetNormal(planeNormal[0], planeNormal[1], planeNormal[2])
plane.SetOrigin(p1[0], p1[1], p1[2])

    # Create the cutter and extract the sampled data
planeCut = vtk.vtkCutter()
planeCut.SetInputData(case.vtkData.VTKObject)
planeCut.GenerateTrianglesOff()
planeCut.SetCutFunction(plane)
planeCut.Update()
cutData = dsa.WrapDataObject(planeCut.GetOutput())