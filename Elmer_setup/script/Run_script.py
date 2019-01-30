#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 14:50:26 2019

@author: Schmulius
"""

import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

ren1 = vtk.vtkRenderer()
ren2 = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren1)
renWin.SetMultiSamples(0)
renWin.AddRenderer(ren2)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)



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

scalar_data = {}

for i in range(narrays):
    
    scalar_data[i] =  vtk_to_numpy(xmlReader.GetOutput().GetPointData().GetArray(i))
    
    print(scalar_data)
        











plane = vtk.vtkPlane()
plane.SetOrigin( 0,0,0 )

volume = vtk.vtkVolume()
volume.SetOrigin( 0,0,0 )
#plane.SetNormal( 0,0,1 )

# pipeline for cutter producing triangles
triCutter = vtk.vtkCutter()
triCutter.SetInputConnection( xmlReader.GetOutputPort() )
triCutter.SetCutFunction( plane )

triMapper = vtk.vtkPolyDataMapper()
triMapper.SetInputConnection( triCutter.GetOutputPort() )
triMapper.ScalarVisibilityOff()

triActor  = vtk.vtkActor()
triActor.SetMapper( triMapper )
triActor.GetProperty().SetColor( 1,0,0 )
triActor.GetProperty().EdgeVisibilityOn()
triActor.GetProperty().SetEdgeColor( 1,1,1 )

#ren1.AddViewProp( triActor )
#ren1.SetViewport( 0,0,0.5,1.0)

# pipeline for cutter producing polygons
polyCutter = vtk.vtkCutter()
polyCutter.GenerateTrianglesOff()
polyCutter.SetInputConnection( xmlReader.GetOutputPort() )
polyCutter.SetCutFunction( plane )

polyMapper = vtk.vtkPolyDataMapper()
polyMapper.SetInputConnection( polyCutter.GetOutputPort() )
polyMapper.ScalarVisibilityOff()

polyActor  = vtk.vtkActor()
polyActor.SetMapper( polyMapper )
polyActor.GetProperty().SetColor( 0,0,1 )
polyActor.GetProperty().EdgeVisibilityOn()
polyActor.GetProperty().SetEdgeColor( 1,1,1 )

ren2.AddViewProp( polyActor )
ren2.SetViewport( 0.5,0,1.0,1.0 )

# the render window
renWin.SetSize(600,500)
iren.Initialize()
