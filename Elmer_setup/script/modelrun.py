45# -*- coding: utf-8 -*-

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
from vtk.util.numpy_support import vtk_to_numpy
#from sh import mount



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
import glob
 


    
'''
Check if file server is connected with local machine

'''

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
    
# Check if volume is mounted, if not mount volume
    
    
    
    def __init__(self, 
                 bump_amplitude, 
                 bump_distribution_x,
                 bump_distribution_y,
                 bump_offset,
                 timestep):
        
        
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
        self.timestep = timestep
        #self.grid_refinement = grid_refinement
        
        
        home_directory = '/Volumes/esd01/docs/jloos/data_small/runs_elmerice_'
        sub_mesh_directory = 'Mesh/'

        # change output directory according to mesh refinement:
        grid_refinement = bool(input("Enter case for folder (True==fixed, False==fixed): "))
        
        if grid_refinement==True:
            res_folder = os.path.join(home_directory + 'fixed')
        else:
            res_folder = os.path.join(home_directory + 'fixed')
            
        # create fodername of the run:    
        run_folder = 'Mesh{:}_{:}{:}_{:}'.format(self.bump_amplitude,self.bump_distribution_x,self.bump_distribution_y,self.bump_offset) 
        # path to the directory of the model run:
        self.run_directory = os.path.join(res_folder,run_folder,sub_mesh_directory)  
        
        
        
        
        
   
    
    
    def print_timestep(self,timestep):
        
        """
        print_timestep 
        
        Parameters
        ----------
        timestep:int
        
        
        This method provides a dictionary of all pvtu files, which where
        treated as timesteps (each pvtu file is equal to one timestep)
        """
        
        self.timestep = timestep
            
        self.dic_timesteps = glob.glob(self.run_directory + '*pvtu')
        self.dic_timesteps.sort(key=os.path.getmtime)
        self.f_name = self.dic_timesteps[self.timestep]
        
        print("This run got", len(self.dic_timesteps), "timesteps.")
        print(self.f_name)
       
        self.path_timestep = os.path.join(self.run_directory, self.f_name)
        self.dict_var = {}
        
        
        #Define xmlreader (file-input for VTU-path)
          
    
        
        return self.f_name
        
        
       
        
    
    def scalar(self):
        pass
        '''
        '''
        #self.scalar= scalar
        
        # define reader for pvtu, set timestape
        xmlReader = vtk.vtkXMLPUnstructuredGridReader()
        xmlReader.SetFileName(self.fname)
        xmlReader.Update()
        
        # number of arrays
        narrays = xmlReader.GetOutput().GetPointData().GetNumberOfArrays()
        # number of partitions
        npartition = xmlReader.GetNumberOfPieces()
        output = xmlReader.GetOutput()
        xmlReader.Update()
        
        
        self.dict_scalars = {}

        for i in range(narrays):
    
                self.dict_scalars[i] = xmlReader.GetPointArrayName(i)
        
        print(self.dict_scalars)
        
        return self.dict_scalars
        
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
        
        def reader(self):
            pass
        
        xmlReader = vtk.vtkXMLPUnstructuredGridReader()
        xmlReader.SetFileName(self.path_timestep)
        xmlReader.Update()
        narrays = xmlReader.GetOutput().GetPointData().GetNumberOfArrays()
        output = xmlReader.GetOutput()
        xmlReader.Update()
        npartition = xmlReader.GetNumberOfPieces()
        
        
        dict_var = {}

        for i in range(narrays):
    
            dict_var[i] = xmlReader.GetPointArrayName(i)
        
        self.velocity =  vtk_to_numpy(xmlReader.GetOutput().GetPointData().GetArray(19))
        
        
        
        print(dict_var)
        print(narrays)
        print(npartition)

   
        
            
            





  
        
        
    
      
      
      
        
        
         
         
     
   
           
           
           
                        
    
            
            
    

