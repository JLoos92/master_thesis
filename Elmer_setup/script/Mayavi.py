#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 12:57:13 2018

@author: Schmulius
"""

from mayavi.app import Mayavi

       

class MyClass(Mayavi):

    def run(self):
        script = self.script
        # `self.script` is the MayaVi Script interface (an instance of
        # enthought.mayavi.script.Script) that is created by the
        # base `Mayavi` class.  Here we save a local reference for
        # convenience.

        ## import any Mayavi modules and filters you want (they must be done here!)
        #script = '/volumes/esd01/docs/jloos/data_small/runs_elmerice_refined/Mesh60_200500_0/Bump_0001.pvtu'
        
        from mayavi.sources.vtk_file_reader import VTKFileReader
        from mayavi.sources.api import VTKXMLFileReader
        from mayavi.modules.scalar_cut_plane import ScalarCutPlane
        from pyface.timer.api import Timer
        from mayavi.core.api import Engine
        from mayavi.sources.vtk_file_reader import VTKFileReader
        from mayavi.modules.surface import Surface

        script.new_scene()                      # to create the rendering scene
  
        
        

if __name__ == '__main__':

    mc = MyClass()
    mc.main()