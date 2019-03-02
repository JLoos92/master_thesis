
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:36:40 2019

@author: Schmulius
"""

from modelrun import *
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
import plotly as py
import numpy as np
from tempfile import TemporaryFile
timestep_path_2 = ModelRun(50,500,500,0,1).print_timestep(0)
timestep_path_60 = ModelRun(50,500,500,0,1).print_timestep(60)




xmlReader2 = vtk.vtkXMLPUnstructuredGridReader()
xmlReader60 = vtk.vtkXMLPUnstructuredGridReader()



xmlReader2.SetFileName(timestep_path_2)
xmlReader60.SetFileName(timestep_path_60)
xmlReader2.Update()
xmlReader60.Update()



npoints=xmlReader2.GetOutput().GetNumberOfPoints()
ncells=xmlReader2.GetOutput().GetNumberOfCells()
#ntuples = xmlReader2.GetOutput().GetNumberOfTuples()





Cells2 =  vtk_to_numpy(xmlReader2.GetOutput().GetCells().GetData())
Cells60 =  vtk_to_numpy(xmlReader60.GetOutput().GetCells().GetData())
xmlReader2.Update()
xmlReader60.Update()

narrays = xmlReader2.GetOutput().GetPointData().GetNumberOfArrays()

dict_var = {}

for i in range(narrays):
    
        dict_var[i] = xmlReader2.GetPointArrayName(i)
        
print(dict_var)

#stress_field2 =  vtk_to_numpy(xmlReader2.GetOutput().GetPointData().GetArray(20))
#stress_field60 =  vtk_to_numpy(xmlReader60.GetOutput().GetPointData().GetArray(20))








