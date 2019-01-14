#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 12:00:05 2018

@author: Schmulius
"""

# Standard library imports
from mayavi.mlab import *
from os.path import join, dirname
import numpy

import sys
import vtk
from mayavi.sources.vtk_file_reader import VTKFileReader
from mayavi.sources.api import VTKXMLFileReader
import os
from pyface.timer.api import Timer
from mayavi.core.api import Engine
from mayavi.sources.vtk_file_reader import VTKFileReader
from mayavi.modules.surface import Surface
from mayavi.modules.scalar_cut_plane import ScalarCutPlane



engine = Engine()
engine.start()
scene = engine.new_scene()


#vtkFile = 'simple.vtk'
vtkFile = '/volumes/esd01/docs/jloos/data_small/runs_elmerice_refined/Mesh60_200500_0/Bump_0001.pvtu'

#reader = VTKFileReader()
#reader = VTKXMLFileReader()

#reader.initialize(vtkFile)
#engine.add_source(reader)

# Add Surface Module
#surface = Surface()
#engine.add_module(surface), 


scp = ScalarCutPlane() # set scp as ScalarCutPlane() module
engine.add_module(scp) # add module to the scene
scp.implicit_plane.normal = (1, 0, 0) # set normal to Ox axis
# set origin to (i=10, j=25, k=25) i.e. integers for a structured grid
scp.implicit_plane.origin = (10, 25, 25)
# set origin to (x=1.0, y=2.5, z=2.5) i.e. reals for unstructured grids
# scp.implicit_plane.origin = (1.0, 2.5, 2.5)
scp.implicit_plane.widget.enabled = False
scp.actor.property.diffuse = 0.0 # set some color properties
scp.actor.property.ambient = 1.0 # 
scp.actor.property.opacity = 1.0 #
scp.module_manager.scalar_lut_manager.data_range = [0, 1]
