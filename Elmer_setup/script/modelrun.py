# -*- coding: utf-8 -*-

"""
Created on Fri Dec 21 12:57:13 2018

@author: Schmulius

This class holds this main properties of the modelpath
"""
import os
import math
import sys
import numpy
import vtk
import subprocess
#from sh import mount



# All returned arrays are cast into either numpy or numarray arrays
#arr=numpy.array








#from timestep import Timestep

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
from mayavi.modules.surface import Surface


# Check if volume is mounted, if not mount volume

path_to_volume = '/Volumes/esd01'

if os.path.ismount(path_to_volume)==True:
        print ('ESD01 fileserver is already mounted.')
else:
        subprocess.Popen(["open smb://134.2.5.43/esd01"], shell=True)
        print ('Fileserver has been mounted.')
        





class ModelRun():
    """
    Class ModelRun
    Here description
    """
    
    def __init__(self, 
                 bump_amplitude, 
                 bump_distribution_x,
                 bump_distribution_y,
                 bump_offset,
                 grid_refinement):
        
        
        """
        ModelRun 
        
        Parameters
        ----------
        bump_amplitude : int
            height of the bump
        bump_distribution_x : int
            height of the bump
        bump_distribution_y : int
            height of the bump    
        bump_offset : int
            extent flow direction
         grid_refinement : boolean
            whether the grid is refined or not
        """
        
        self.bump_amplitude = bump_amplitude
        self.bump_distribution_x = bump_distribution_x
        self.bump_distribution_y = bump_distribution_y
        self.bump_offset = bump_offset
        self.grid_refinement = grid_refinement
        
        
        home_directory = '/volumes/esd01/jloos/data_small/runs_elmerice_'
        sub_mesh_directory = 'Mesh/'

        # change output directory according to mesh refinement:
        grid_refinement = bool(input("Enter case for folder (True==refined, False==uniform): "))
        
        if grid_refinement==True:
            res_folder = os.path.join(home_directory + 'refined')
        else:
            res_folder = os.path.join(home_directory + 'uniform')
            
        # create fodername of the run:    
        run_folder = 'Mesh{:}_{:}{:}_{:}'.format(bump_amplitude,bump_distribution_x,bump_distribution_y,bump_offset) 
        # path to the directory of the model run:
        self.run_directory = os.path.join(res_folder,run_folder,sub_mesh_directory)       
        print(self.run_directory)
        
       # model_path = self.run_directory

        
        
     
     
           
           
           
           
                        
    
            
            
    

