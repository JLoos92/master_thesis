#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 12:57:13 2018

@author: Schmulius
"""

#__________________________________ modules & packages ____________________________________#

from mayavi.modules.scalar_cut_plane import ScalarCutPlane
from pyface.timer.api import Timer
from mayavi.core.api import Engine
from mayavi.core.base import Base
from mayavi.core.module import Module
from mayavi.core.lut_manager import LUTManager
from mayavi.core.common import handle_children_state, exception
from mayavi.core.pipeline_info import PipelineInfo
from mayavi import mlab
from mayavi import tools     
from numpy import array


#---------------------------------------------------------------------------------------
#__________________________________ start mayavi engine _______________________________#
#---------------------------------------------------------------------------------------
try:
    engine = mayavi.engine
except NameError:
    from mayavi.api import Engine
    engine = Engine()
    engine.start()
if len(engine.scenes) == 0:
    engine.new_scene()
    
    
    
#---------------------------------------------------------------------------------------    
#__________________________________ load data ____________________________________#
#---------------------------------------------------------------------------------------
    
vtkxml_file_reader = engine.open('/Volumes/esd01/docs/jloos/data_small/runs_elmerice_refined/Mesh20_5001000_0_0.1/Mesh/ForwardRemesh100M20L0011.pvtu')
# Set time-step for visualisation; depending on simulation setting time-step cann range between 
# 0.1 and 0.5 = 2 (two time sizes 1a)
from mayavi.modules.surface import Surface
surface = Surface()
engine.add_module(surface, obj=None)



# Set timestep for visualisation
# only depends on pvtu, choose random pvtu and then set time step
#vtkxml_file_reader = engine.scenes[0].children[0]
#vtkxml_file_reader.name = 'VTK XML file (ForwardRemesh100M20L0004.pvtu) (timeseries)'
#vtkxml_file_reader.file_path = '/Volumes/esd01/docs/jloos/data_small/runs_elmerice_refined/Mesh20_5001000_0_0.1/Mesh/ForwardRemesh100M20L0004.pvtu'
vtkxml_file_reader.timestep = 3




#---------------------------------------------------------------------------------------
#_______________________________ Set scene_______________________________#
#---------------------------------------------------------------------------------------

# Arrange dataset and frame for scene. Choose between set of point scalar names
# (ds/dt, groundedmask,...)and point vector names (velocity,,...)
# set scale(aspect ratio) between 20-30 for better visualisation


vtkxml_file_reader.point_scalars_name = 'groundedmask'


#surface = engine.scenes[0].children[0].children[0].children[0]
surface.actor.actor.orientation = array([ 0., -0.,  0.])
surface.actor.actor.origin = array([0., 0., 0.])
surface.actor.actor.position = array([0., 0., 0.])
surface.actor.actor.scale = array([ 1.,  1., 20.])
surface.actor.actor.render_time_multiplier = 0.6115716831980361
surface.actor.actor.scale = array([ 1.,  1., 20.])






import numpy
from mayavi.mlab import *

def test_volume_slice():
    x, y, z = np.ogrid[-5:5:64j, -5:5:64j, -5:5:64j]

    scalars = x * x * 0.5 + y * y + z * z * 2.0

    obj = volume_slice(scalars, plane_orientation='x_axes')
    return obj






