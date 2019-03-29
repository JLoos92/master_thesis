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
import glob
from numpy import array
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
import concave_hull
import numpy as np
from tempfile import TemporaryFile


#from timestep import Timestep

#from mayavi.modules.scalar_cut_plane import ScalarCutPlane
#from pyface.timer.api import Timer
#from mayavi.core.api import Engine
#from mayavi.core.base import Base
#from mayavi.core.module import Module
#from mayavi.core.lut_manager import LUTManager
#from mayavi.core.common import handle_children_state, exception
#from mayavi.core.pipeline_info import PipelineInfo
#from mayavi import mlab
#from mayavi import tools     

#from mayavi.modules.surface import Surface
import glob
 


    
'''
Check if file server is connected with local 
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
       # grid_refinement = bool(input("Enter case for folder (True==fixed, False==fixed): "))
        
        #if grid_refinement==True:
         #   res_folder = os.path.join(home_directory + 'fixed')
        #else:
        res_folder = os.path.join(home_directory + 'fixed')
            
        # create fodername of the run:    
        run_folder = 'Mesh{:}_{:}{:}_{:}'.format(self.bump_amplitude,self.bump_distribution_x,self.bump_distribution_y,self.bump_offset) 
        # path to the directory of the model run:
        self.run_directory = os.path.join(res_folder,run_folder,sub_mesh_directory)  
        
        
        
        """        
        This method provides a dictionary of all pvtu files, which where
        treated as timesteps (each pvtu file is equal to one timestep)
        """
        
        self.timestep = timestep
        self.dic_timesteps = glob.glob(self.run_directory + '*pvtu')
        self.dic_timesteps.sort(key=os.path.getmtime)
        self.f_name = self.dic_timesteps[self.timestep]
        
        print("Timesteps = ", len(self.dic_timesteps))
        print(self.f_name)
       
        self.path_timestep = os.path.join(self.run_directory, self.f_name)
        self.dict_var = {}
        
        
        # Define xmlreader (file-input for VTU-path)
        # Create vtu objects and to collect properties
        
        self.xmlReader = vtk.vtkXMLPUnstructuredGridReader()
        self.xmlReader.SetFileName(self.f_name)
        self.xmlReader.Update()
        
        # Get outputs
        self.npoints = self.xmlReader.GetOutput().GetNumberOfPoints()
        self.ncells = self.xmlReader.GetOutput().GetNumberOfCells()
        self.bounds = self.xmlReader.GetOutput().GetBounds()
        self.narrays = self.xmlReader.GetOutput().GetPointData().GetNumberOfArrays()
        self.Points = vtk_to_numpy(self.xmlReader.GetOutput().GetPoints().GetData())
        self.Cells =  vtk_to_numpy(self.xmlReader.GetOutput().GetCells().GetData())
        
        
        # Dictionary for all arrays, scalars etc.
        self.dict_var = {}
        #self.dict_arr = {}
        
        for i in range(self.narrays):
    
#            self.dict_var[i] = self.xmlReader.GetPointArrayName(i)
#            self.dict_arr[i] = self.xmlReader.GetOutput().GetPointData().GetArray(i))
            self.dict_var[self.xmlReader.GetPointArrayName(i)] = self.xmlReader.GetOutput().GetPointData().GetArray(i)
        # Pick needed array
       
        
       
                
        #iteration über 
       # for key in dict_var.keys:
           # dict_var[key]
        #return self.xmlReader

       
        
    
    def cutter(self, GL=None):
        
        """
        Method: cutter
        
        Parameters
        ----------
        -
        
        This method performes cutting for given run (default is cutting plane normal
        in yz with origin at GL. Further methods and functions relate to the original vtk-class 
        or the ModelRun-class).
        Cutting provides also profile-lines for extracting intercepting values of given array.
        
        """
        
        
        
        # Change for cut plane in x direction ()
        if GL is None:
            self.GL = 1065000
        elif GL is not None:
            self.GL = GL

        #create a plane to cut,here it cuts in the XZ direction (xz normal=(1,0,0);XY =(0,0,1),YZ =(0,1,0)
        self.plane=vtk.vtkPlane()
        self.plane.SetOrigin(self.GL,0,0)
        self.plane.SetNormal(1,0,0)

        #create cutter
        self.clipData = vtk.vtkClipDataSet()   
        self.clipData.SetClipFunction(self.plane)
        self.clipData.SetInputConnection(self.xmlReader.GetOutputPort())
        self.clipData.Update()
        #self.clipData = self.clipData.GetOutput()
       # self.dict_var_clipped = self.clipData.SetInputConnection(self.xmlReader.GetOutputPort().GetPointArrayName())
        
        #Check for arrays
       # for i in range(self.narrays):
    
           
           # self.dict_var[self.clipData.GetPointArrayName(i)] = self.clipData.GetOutput().GetPointData().GetArray(i)
        
       
        
        return self.clipData
        

    
        
    def convexhull(self):
        
        """
        Method: convexhull 
        
        Parameters
        ----------
        self
        
       
        
        """
       # self.GetOutput = self.cutter().GetOutput()
        self.triangulation = vtk.vtkDelaunay3D()
        self.triangulation.SetInputData(self.cutter().GetOutput())
        self.triangulation.Update()
        
        self.convexhull = vtk.vtkDataSetSurfaceFilter() 
        self.convexhull.SetInputConnection(self.triangulation.GetOutputPort())
        self.convexhull.Update()
        
        self.ch_points = self.convexhull.GetOutput()
        
        return self.ch_points
              

            

    def compute_hydrostatic_thickness(self):
        
        """
        
        Method: compute_hydrostatic_thickness 
        ----------
        This method computes the hydrostatic thickness for a clipped area,
        which is defined by the cutter method.
        
        Parameters
        ----------
        -
        
        Returns:
        ----------
        Tuple(float): x-coordinates, y-coordinates, hydrostatic thickness
        
        
        """
     
    # Paraneter for calculation of hydrostatic thickness
        
        p_w = 1027.0 #kg m−3 ), ice (ρi =918kgm−3), and air (ρa =2kgm−3):
        p_i = 900.0
        p_a = 2.0
        H_a = 0
        
        #self.cutter = self.clipData
        self.clipped_area = self.cutter()
        self.zs = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray('zs'))
        self.zb = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray('zb'))
        self.points = vtk_to_numpy(self.clipped_area.GetOutput().GetPoints().GetData())
        
        
       
        self.x = self.points[:,0]
        self.y = self.points[:,1]
        self.z = self.points[:,2]

       
        def selectMinz_vectorized(x, y, z):
            # Get grouped lex-sort indices
            sidx = (y + x*(y.max() - y.min() + 1)).argsort()
            # or sidx = np.lexsort([x, y])
        
            # Lex-sort x, y, z
            x_sorted = x[sidx]
            y_sorted = y[sidx]
            z_sorted = z[sidx]
        
            # Get equality mask between each sorted X and Y elem against previous ones.
            # The non-zero indices of its inverted mask gives us the indices where the 
            # new groupings start. We are calling those as cut_idx.
            seq_eq_mask = (x_sorted[1:] == x_sorted[:-1]) & (y_sorted[1:] == y_sorted[:-1])
            cut_idx = np.flatnonzero(np.concatenate(( [True], ~seq_eq_mask)))
        
            # Use those cut_idx to get intervalled minimum values
            minZ = np.minimum.reduceat(z_sorted, cut_idx)
        
            # Make tuples of the groupings of x,y and the corresponding min Z values
            return (zip(x_sorted[cut_idx], y_sorted[cut_idx]), minZ.tolist())
        
        
        def selectMaxz_vectorized(x, y, z):
            # Get grouped lex-sort indices
            sidx = (y + x*(y.max() - y.min() + 1)).argsort()
            # or sidx = np.lexsort([x, y])
        
            # Lex-sort x, y, z
            x_sorted = x[sidx]
            y_sorted = y[sidx]
            z_sorted = z[sidx]
        
            # Get equality mask between each sorted X and Y elem against previous ones.
            # The non-zero indices of its inverted mask gives us the indices where the 
            # new groupings start. We are calling those as cut_idx.
            seq_eq_mask = (x_sorted[1:] == x_sorted[:-1]) & (y_sorted[1:] == y_sorted[:-1])
            cut_idx = np.flatnonzero(np.concatenate(( [True], ~seq_eq_mask)))
        
            # Use those cut_idx to get intervalled maximum values
            maxZ = np.maximum.reduceat(z_sorted, cut_idx)
        
            # Make tuples of the groupings of x,y and the corresponding min Z values
            return (zip(x_sorted[cut_idx], y_sorted[cut_idx]), maxZ.tolist())
        
        
        zmin = selectMinz_vectorized(self.x,self.y,self.z)
        zmax = selectMaxz_vectorized(self.x,self.y,self.z)

        zmax = zmax[1]
        zmax = np.asarray(zmax)

        zmin = zmin[1]
        zmin = np.asarray(zmin)
        
        
        zmin_ind = np.isin(self.z,zmin)
        zmax_ind = np.isin(self.z,zmax)
        
        zmin_ind = np.where(zmin_ind)
        zmax_ind = np.where(zmax_ind)
        
        self.zb_new = self.zb[zmin_ind]
        self.zs_new = self.zs[zmax_ind]

#        self.zs_ind = np.asarray(np.where(self.zs>0))
#        self.zs_ind = self.zs_ind.transpose().flatten()
        self.x_new = self.x[zmax_ind]
        self.y_new = self.y[zmax_ind]
#        
#        self.zs_s = self.zs[self.zs!=0]
#        self.zb_s = self.zb[self.zb!=0]
        
        self.thick_model = -self.zs_new-self.zb_new       
        self.thick_calc = np.divide((p_w*self.zs_new),(p_w-p_i))                
        self.h_thickness = self.thick_model - self.thick_calc 
        
        return self.x_new,self.y_new,self.h_thickness
    
    
    
    
    
    
    def cut_and_slice(self, 
                      xcoord,
                      var):
        
        """
        
        Method: cut_and_slice
        ----------
        This methods cuts the given domain at a x-coordinate, which is defined
        by the method itself. The Variable input defines the output array and
        is defined and listed in dict_var.
        (Either scalar or vector)
        
        Parameters
        ----------
        xcoord : int
        var : 'string' (execute ModelRun with dict_var to  check for suitable input)
        
        Returns:
        ----------
        Tuple(float): Matrix of variable cut at x-coordinate, line points to
        compute concave hull
        
        
        """
        # Set xcoord and variable for cut and slice method 
        
        self.xcoord = xcoord
        self.var = var
        
        #Cutter1
        
        self.line=vtk.vtkPlane()
        self.line.SetOrigin(self.xcoord,0,0)
        self.line.SetNormal(1,0,0)
        
        self.clipped_area = self.cutter()
        self.cutter = vtk.vtkCutter()   
        self.cutter.SetCutFunction(self.line)
        self.cutter.SetInputConnection(self.clipped_area.GetOutputPort())
        self.cutter.Update()
        self.line_points = vtk_to_numpy(self.cutter.GetOutput().GetPoints().GetData())
        
        self.out = vtk_to_numpy(self.cutter.GetOutput().GetPointData().GetArray(self.var))


        return self.out, self.line_points


        
 
    def compute_concavehull(self,
                            xcoord,
                            neighbors):
       # pass
        
        """
        Method: compute_concavehull 
        
        Parameters
        ----------

        
       
        
        """     
        # Set x-coordinate
        self.xcoord = xcoord
        
        # Set number for neighbors to compute concave hull
        self.neighbors = neighbors
        
        
        # Cut and slice
       
        self.line=vtk.vtkPlane()
        self.line.SetOrigin(self.xcoord,0,0)
        self.line.SetNormal(1,0,0)
        
        self.clipped_area = self.cutter()
        self.cutter = vtk.vtkCutter()   
        self.cutter.SetCutFunction(self.line)
        self.cutter.SetInputConnection(self.clipped_area.GetOutputPort())
        self.cutter.Update()
        self.line_points = vtk_to_numpy(self.cutter.GetOutput().GetPoints().GetData())
        
        # Rearrange matrix to fit input shape ()
        self.hull_points = np.delete(self.line_points, 0, 1)
        
        # Compute concave hull at crosssection of line points
        self.hull = concave_hull.compute(self.hull_points,self.neighbors)
        
        return self.hull
  
        
        
    
      
      
      
        
        
         
         
     
   
           
           
           
                        
    
            
            
    

