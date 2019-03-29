#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 12:59:31 2018

@author: Schmulius
"""
from mayavi.modules.scalar_cut_plane import ScalarCutPlane

from pyface.timer.api import Timer
from mayavi.core.api import Engine

from mayavi.core.base import Base
from mayavi.core.module import Module
from mayavi.core.lut_manager import LUTManager
from mayavi.core.common import handle_children_state, exception
from mayavi.core.pipeline_info import PipelineInfo



scp = Engine()
scp.start()
scene = scp.new_scene()


scp = ScalarCutPlane() # set scp as ScalarCutPlane() module
 # add module to the scene
scp.implicit_plane.normal = (1, 0, 0) # set normal to Ox axis
# set origin to (i=10, j=25, k=25) i.e. integers for a structured grid
scp.implicit_plane.origin = (10, 25, 25)
# set origin to (x=1.0, y=2.5, z=2.5) i.e. reals for unstructured grids
# scp.implicit_plane.origin = (1.0, 2.5, 2.5)
scp.implicit_plane.widget.enabled = True
scp.actor.property.diffuse = 0.0 # set some color properties
scp.actor.property.ambient = 1.0 # 
scp.actor.property.opacity = 1.0 #
#scp.module_manager.scalar_lut_manager.data_range = [0, 1]
        