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





#Section  with VTK-Files




reader = vtk.vtkDataSetReader()
reader.SetFileName("/Volumes/esd01/docs/datasmall/runselmerice/Mesh200_5001000_0/Forward5KM10L0001par0001.vtu")
reader.ReadAllScalarsOn()  # Activate the reading of all scalars
reader.Update()
data=reader.GetOutput()