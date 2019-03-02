#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 14:50:26 2019

@author: Schmulius
"""

import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy




# Setup Reader and path
xmlReader = vtk.vtkXMLPUnstructuredGridReader()
xmlReader.SetFileName('/Volumes/esd01/docs/jloos/data_small/runs_elmerice_refined/Mesh20_5001000_0_0.1/Mesh/ForwardRemesh100M20L0011.pvtu')
xmlReader.Update()
narrays = xmlReader.GetOutput().GetPointData().GetNumberOfArrays()
output = xmlReader.GetOutput()
xmlReader.Update()
npartition = xmlReader.GetNumberOfPieces()
#bedrock = xmlReader.GetPointArrayName(4)

# Create variable-dictionary (contains scalars and vector)

dict_var = {}

for i in range(narrays):
    
        dict_var[i] = xmlReader.GetPointArrayName(i)
        
print(dict_var)
    
xmlReader.Update()
npoints=xmlReader.GetOutput().GetNumberOfPoints()
ncells=xmlReader.GetOutput().GetNumberOfCells()


xmlReader.Update()
Cells =  vtk_to_numpy(xmlReader.GetOutput().GetCells().GetData())
xmlReader.Update()
Points = vtk_to_numpy(xmlReader.GetOutput().GetPoints().GetData())


    
velocity =  vtk_to_numpy(xmlReader.GetOutput().GetPointData().GetArray(18))

test = xmlReader.GetOutput().GetPointData().SetActiveScalars("velocity")











  












