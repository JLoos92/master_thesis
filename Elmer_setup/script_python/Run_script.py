#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 14:50:26 2019

@author: Schmulius
"""

import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy




p1 = (1010000,0)
p2 = (1080000,0)

p1 = np.append(np.array(p1),0)
p2 = np.append(np.array(p2),0)


unit = (p2-p1/np.linalg.norm(p2-p1))

geoNormal = np.array([0,0,1])
planeNormal = np.cross(unit,geoNormal)

plane = vtk.vtkPlane()
plane.SetNormal(planeNormal[0],planeNormal[0])




  












