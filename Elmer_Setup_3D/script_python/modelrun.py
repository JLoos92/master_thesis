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
import glob
from scipy.interpolate import interp1d
import pandas as pd



 


    
#==============================================================================    
# Check if volume is mounted, if not mount volume
#==============================================================================


path_to_volume = '/Volumes/esd01'

if os.path.ismount(path_to_volume)==True:
        print ('ESD01 fileserver is already mounted.')
else:
        subprocess.Popen(["open smb://134.2.5.43/esd01"], shell=True)
        print ('Fileserver has been mounted to local machine.')  






class ModelRun():
    
    """
    Class ModelRun:
    The ModelRun class simplifies the data arrangement for large, three
    dimensional, parallel .vtu files. This class contains static methods and 
    user functions to load data files of computed model runs of the esd 
    fileserver (or path to be defined where simulation files are saved).
    """

    
    
    def __init__(self, 
                 bump_amplitude, 
                 bump_distribution_x,
                 bump_distribution_y,
                 prop,
                 timestep):
        
        
        """
        Initial: ModelRun 
        
        Parameters
        ----------
        bump_amplitude : int
            height of the bump
        bump_distribution_x : int
            height of the bump
        bump_distribution_y : int
            height of the bump    
        prop : string
            extra specification e.g. double (two gauss functions)
         timestep : int
            number of timestep
        """
        
        self.bump_amplitude = bump_amplitude
        self.bump_distribution_x = bump_distribution_x
        self.bump_distribution_y = bump_distribution_y

        self.prop = prop # extra property e.g. double bump
        self.timestep = timestep
        
        
        
        home_directory = '/Volumes/esd01/docs/jloos/data_small/runs_elmerice_'
        sub_mesh_directory = 'Mesh/'

        # Change output directory according to mesh refinement:
        # grid_refinement = bool(input("Enter case for folder 
        # (True==fixed, False==fixed): "))
        
        #if grid_refinement==True:
         #   res_folder = os.path.join(home_directory + 'fixed')
        #else:
        
        self.res_folder = os.path.join(home_directory + 'fixed')
            
        # create fodername of the run:    
        run_folder = 'Mesh{:}_{:}{:}_{:}'.format(self.bump_amplitude,
<<<<<<< HEAD:Elmer_Setup_3D/script_python/modelrun.py
                          self.bump_distribution_x,
                          self.bump_distribution_y,
                          self.prop) 
        
        # path to the directory of the model run:
        self.run_directory = os.path.join(self.res_folder,
                                          run_folder,
=======
                          self.bump_distribution_x,self.bump_distribution_y,
                          self.prop) 
        
        # path to the directory of the simulation:
        self.run_directory = os.path.join(self.res_folder,run_folder,
>>>>>>> d0d26ebd5636d5e00ad413393f618c412a6dc2df:Elmer_setup/script_python/modelrun.py
                                          sub_mesh_directory)  
        
        
        
        #======================================================================        
        # The following segment provides a dictionary of all .pvtu files
        # and a subdirectory list of all runs. 
        #====================================================================== 
        
        dirlist = [item for item in os.listdir(self.res_folder) 
                  if os.path.isdir(os.path.join(self.res_folder,item))]
        
        self.dirlist = dirlist
                    
       
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
        self.narrays = self.xmlReader.GetOutput().GetPointData(). \
                                    GetNumberOfArrays()
        self.Points = vtk_to_numpy(self.xmlReader.GetOutput(). \
                                   GetPoints().GetData())
        self.Cells =  vtk_to_numpy(self.xmlReader.GetOutput(). \
                                   GetCells().GetData())
        
        
        # Dictionary for all arrays, scalars etc.
        self.dict_var = {}
       
        
        for i in range(self.narrays):
    
        #self.dict_var[i] = self.xmlReader.GetPointArrayName(i)
        #self.dict_arr[i] = self.xmlReader.GetOutput().GetPointData().GetArray(i))
            self.dict_var[self.xmlReader.GetPointArrayName(i)] \
            = self.xmlReader.GetOutput().GetPointData().GetArray(i)
        

    
    def cutter(self ,GL=None):
        
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
            self.GL = 1060000
        elif GL is not None:
            self.GL = GL
        
        
        
        # create a plane to cut; here it cuts in the XZ direction
        # (xz normal=(1,0,0);XY =(0,0,1),YZ =(0,1,0)
        self.plane_shelf =vtk.vtkPlane()
        self.plane_shelf.SetOrigin(self.GL,0,0)
        self.plane_shelf.SetNormal(1,0,0)

        # create first cutter
        self.clipData_shelf = vtk.vtkClipDataSet()   
        self.clipData_shelf.SetClipFunction(self.plane_shelf)
        self.clipData_shelf.SetInputConnection(self.xmlReader.GetOutputPort())
        self.clipData_shelf.Update()
        
        # create plane to cut; end of shelf is cut due to boundary and mesh-
        # refinement adjustment for plotting
        self.plane_shelf_end = vtk.vtkPlane()
        self.plane_shelf_end.SetOrigin(1079000,0,0)
        self.plane_shelf_end.SetNormal(-1,0,0)
        
        # create second cutter
        self.clipData = vtk.vtkClipDataSet()   
        self.clipData.SetClipFunction(self.plane_shelf_end)
        self.clipData.SetInputConnection(self.clipData_shelf.GetOutputPort())
        self.clipData.Update()
        
        
        #self.clipData = self.clipData.GetOutput()
       # self.dict_var_clipped = self.clipData.SetInputConnection(self.xmlReader.GetOutputPort().GetPointArrayName())
        #self.out = vtk_to_numpy(self.cutter.GetOutput().GetPointData().GetArray(self.var))
        
        #Check for arrays
        #for i in range(self.narrays):
          
           # self.dict_var[self.clipData.GetPointArrayName(i)] = self.clipData.GetOutput().GetPointData().GetArray(i)
        
       
        
        return self.clipData
        

    
        
    def compute_convexhull(self):
        
        """
        Method: convexhull 
        ----------
<<<<<<< HEAD:Elmer_Setup_3D/script_python/modelrun.py
        Computes boundary hull of cutted domain. Only valid if boundaries
        are rectangular (e.g. hull of 3D- Model)
=======
        Computes boundary hull of cutted domain. This static method can be used 
        to obtain the rectangular boundaries of a three dimensional domain.
>>>>>>> d0d26ebd5636d5e00ad413393f618c412a6dc2df:Elmer_setup/script_python/modelrun.py
        
        Parameters
        ----------
        -
<<<<<<< HEAD:Elmer_Setup_3D/script_python/modelrun.py
        
        Returns:
        ----------
        Tuple(float): Points of convex-hull (x, y, z)
                      
        """
        
        
=======
               
        """
       
        # make triangulation for cutted domain
>>>>>>> d0d26ebd5636d5e00ad413393f618c412a6dc2df:Elmer_setup/script_python/modelrun.py
        self.triangulation = vtk.vtkDelaunay3D()
        self.triangulation.SetInputData(self.cutter().GetOutput())
        self.triangulation.Update()
        
        # apply vtk surface filter with triangulation to obtain hull
        self.convexhull = vtk.vtkDataSetSurfaceFilter() 
        self.convexhull.SetInputConnection(self.triangulation.GetOutputPort())
        self.convexhull.Update()
        
        self.ch_points = self.convexhull.GetOutput()
        
        return self.ch_points
              

 
                  
    def selectMinz(self, x, y, z):
<<<<<<< HEAD:Elmer_Setup_3D/script_python/modelrun.py
        
        """
        Method: selectMinz
        ----------
        Calculates bottom surface by vectorization of 3D- point cloud (e.g. zb)
        
        Parameters
        ----------
        x : numpy-array
            Coordinate of cutted domain    
        
        y : numpy-array
            Coordinate of cutted domain    
        
        z : numpy-array
            Coordinate of cutted domain
        
        Returns:
        ----------
        Tuple(float): x(sorted), y (sorted), min z (sorted)
                
        """
        
=======
        # Get grouped  indices
>>>>>>> d0d26ebd5636d5e00ad413393f618c412a6dc2df:Elmer_setup/script_python/modelrun.py
        sidx = (y + x*(y.max() - y.min() + 1)).argsort()
        
    
        # Sort x, y, z
        x_sorted = x[sidx]
        y_sorted = y[sidx]
        z_sorted = z[sidx]
    
        # Get equality mask between each sorted X and Y elem against previous 
        seq_eq_mask = (x_sorted[1:] == x_sorted[:-1]) & (y_sorted[1:] == y_sorted[:-1])
        cut_idx = np.flatnonzero(np.concatenate(( [True], ~seq_eq_mask)))
    
        # Use those cut_idx to get intervalled minimum values
        minZ = np.minimum.reduceat(z_sorted, cut_idx)
    
        # Make tuples of the groupings of x,y and the corresponding min Z values
        return (x_sorted[cut_idx], y_sorted[cut_idx]), minZ.tolist()
    
         # Function to find min-values of z_s
         
   
    
    def selectMaxz(self,x, y, z):
        
        """
        Method: selectMinz
        ----------
        Calculates top surface by vectorization of 3D - point cloud (eg. zs)
        
        Parameters
        ----------
        x : numpy-array
            Coordinate of cutted domain    
        
        y : numpy-array
            Coordinate of cutted domain    
        
        z : numpy-array
            Coordinate of cutted domain
        
        Returns:
        ----------
        Tuple(float): x(sorted), y (sorted), max z (sorted)
        
        """
        
        sidx = (y + x*(y.max() - y.min() + 1)).argsort()
       
    
        # sort x, y, z
        x_sorted = x[sidx]
        y_sorted = y[sidx]
        z_sorted = z[sidx]
    
        # Get equality mask between each sorted X and Y elem against previous 
        seq_eq_mask = (x_sorted[1:] == x_sorted[:-1]) & (y_sorted[1:] == y_sorted[:-1])
        cut_idx = np.flatnonzero(np.concatenate(( [True], ~seq_eq_mask)))
    
        # Use those cut_idx to get intervalled maximum values
        maxZ = np.maximum.reduceat(z_sorted, cut_idx)
    
        # Make tuples of the groupings of x,y and the corresponding min Z values
        return (x_sorted[cut_idx], y_sorted[cut_idx]), maxZ.tolist()
    
 


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
     
        # Paraneter setup for calculation of hydrostatic thickness
        # !!!!!!!!!! Must be changed if input.sif file is changed!!!!!!!!!!!!!!
        self.p_w = 1000.0 # kg m−3 )
        self.p_i = 900.0  # ice (ρi =918kgm−3)
        p_a = 2.0         # air (ρa =2kgm−3)
        H_a = 0
        
        #self.cutter = self.clipData
        self.clipped_area = self.cutter()
        self.zs = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray('zs'))
        self.zb = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray('zb'))
        self.points = vtk_to_numpy(self.clipped_area.GetOutput().GetPoints().GetData())
        
        
       
        self.x = self.points[:,0]
        self.y = self.points[:,1]
        self.z = self.points[:,2]

        # Function to find min-values of z_b       
        zmin = self.selectMinz(self.x,self.y,self.z)
        zmax = self.selectMaxz(self.x,self.y,self.z)

        zs = zmax[1]
        zs = np.asarray(zs)
        xy_zmax = zmax[0]
       
        
        zb = zmin[1]
        zb = np.asarray(zb)
        xy_zmin = zmin[0]
        

        self.zb_new = zb
        self.zs_new = zs
        
        self.x_corr = xy_zmax[0] 
        self.y_corr = xy_zmax[1]
        
        # Calculation of real hydrostatic thickniss, modelled thickniss and
        # deviation of hydrostatic equilibrium
        self.thick_model = -self.zs_new+self.zb_new       
        self.thick_calc = np.divide((self.p_w*self.zs_new),(self.p_w-self.p_i))                
        self.h_thickness = self.thick_model + self.thick_calc 
        
        
        
        # Return values of coordinates and calculated hydrostatic thickness
        return self.x_corr, self.y_corr, self.h_thickness, self.thick_model
    
    
    
    
    
    
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
        Tuple (float): Matrix of variable cut at x-coordinate, line points to
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
                            xcoord):
       # pass
        
        """
        Method: compute_concavehull 
        ----------
        This method computes the channel imprint at a given x - coordinate.
        Define number of neighbors (default = 3) for concave hull function and
        x-coordinate.
        
        Parameters
        ----------
        xcoord : int
        
        
        Returns:
        ----------
        Tuple (float): x, z arrays for concave hull
        
        """     
        
        
        # Set x-coordinate
        self.xcoord = xcoord

        
        
        # Cut and slice
        # for x in self.coordsx:
            
        self.line=vtk.vtkPlane()
        self.line.SetOrigin(xcoord,0,0)
        self.line.SetNormal(1,0,0)
        
        self.clipped_area = self.cutter()
        self.cutter_ch = vtk.vtkCutter()   
        self.cutter_ch.SetCutFunction(self.line)
        self.cutter_ch.SetInputConnection(self.clipped_area.GetOutputPort())
        self.cutter_ch.Update()
        self.line_points = vtk_to_numpy(self.cutter_ch.GetOutput().GetPoints().GetData())

       
        # Rearrange matrix to fit input shape ()
        self.hull_points = np.delete(self.line_points, 0, 1)
        
       
        
        #======================================================================
        # Crosssection for zs ()
        #======================================================================
        
        #Cutter for hydrostatic thickness section        
        self.line_h=vtk.vtkPlane()
        self.line_h.SetOrigin(self.xcoord,0,0)
        self.line_h.SetNormal(1,0,0)
        
        self.clipped_area = self.cutter()
        self.cutter_line = vtk.vtkCutter()   
        self.cutter_line.SetCutFunction(self.line_h)
        self.cutter_line.SetInputConnection(self.clipped_area.GetOutputPort())
        self.cutter_line.Update()
        self.line_points_calc = vtk_to_numpy(self.cutter_line.GetOutput().GetPoints().GetData())
        
        self.zmax = self.selectMaxz(self.line_points_calc[:,0],\
                          self.line_points_calc[:,1],self.line_points_calc[:,2])
       
        
        xy_zmax = self.zmax[0]
       
        x  = xy_zmax[0] 
        y  = xy_zmax[1]
        
        
                
        zs_line = self.zmax[1]
        zs_line = np.asarray(zs_line)
        
        
        self.thick_calc = np.divide((self.p_w*zs_line),(self.p_w-self.p_i))
        self.thick_calc = self.thick_calc *(-1)
        self.thick_calc = self.thick_calc + zs_line
        
        self.hull_points_calc =  y, self.thick_calc
        self.hull_points_calc = np.asarray(self.hull_points_calc).T
        
        
               
        
        # Using pandas-dataframe
        df = pd.DataFrame(self.hull_points,columns=['x','y'])
        df_calc = pd.DataFrame(self.hull_points_calc,columns=['x','y'])
                      
        # Group
        grouped = df.groupby('x')
        grouped_calc = df_calc.groupby('x')
        
        #Upper and lower boundary
        lower_boundary = grouped.min()
        upper_boundary = grouped.max()
        
        #Lower boundary for calculated hydrostatic thickness
        original_thickness = grouped_calc.min()
        
        
        self.hull_panda = upper_boundary,lower_boundary,original_thickness
        
        # Add function for smoothing concave hull
            
               
        
        return upper_boundary, lower_boundary, original_thickness
  
        
        
    
      
      
      
        
        
         
         
     
   
           
           
           
                        
    
            
            
    

