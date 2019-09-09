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

import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'


 


    
#==============================================================================    
# ESD01 fileserver
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
    
    """

    
    
    def __init__(self, 
             bump_amplitude , 
             bump_distribution_x ,
             bump_distribution_y ,
             prop,
             timestep ,
             dimensions = None,
             **kwargs):
        
        
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
            extra specification e.g. double (two gauss functions in 3d)
            or grid extent in 2d
        
        timestep : int
            number of timestep
        
        dimensions : string
            2-D or 3-D
        """
        
        
        
        
        
        # define paths to destination of output       
        self.home_directory = '/Volumes/esd01/docs/jloos/data_small/runs_elmerice_'
        self.sub_mesh_directory = 'Mesh/'
        self.sub_mesh_directory_2d = 'channel2d/'
        
        
        
        # change output directory according to dimensions (2d or 3d):
        self.bump_amplitude = bump_amplitude
        self.bump_distribution_x = bump_distribution_x
        self.bump_distribution_y = bump_distribution_y
        self.prop = prop # extra property e.g. double bump
        self.timestep = timestep
        self.dimensions = dimensions
        
        if self.dimensions is None:
             self.res_folder = os.path.join(self.home_directory + 'fixed')
             print('Made ModelRun object for 3D-case')
             
             
        elif self.dimensions==str("2"):            
             self.res_folder = os.path.join(self.home_directory + '2d')
             print('Made ModelRun object for 2D-case.')
        else:
            raise ValueError('This filetype does not exist in the simulation'
            'folder. Check if dimension input parameter is set to 2 or to None for 2d or 3d case ,'
            'respectively.')
        


        
        #======================================================================  
        # Paraneter setup for calculation of hydrostatic thickness
        # !!!!!!!! Must be changed if input.sif file is changed!!!!!!!!
        #====================================================================== 
        
        self.p_w = 1000.0 # kg m−3 )
        self.p_i = 900.0  # ice 
        self.p_a = 2.0    # air (ρa =2kgm−3)
        self.H_a = 0      # firn corrections              
                  
        #======================================================================        
        # The following segment provides a dictionary of all .pvtu or vtu files
        # (depending on 2d or 3d case) and a subdirectory list of all runs. 
        #====================================================================== 
        
        
        # Dictionary of folder-names (runs) for
        dirlist = [item for item in os.listdir(self.res_folder) \
                  if os.path.isdir(os.path.join(self.res_folder,item))]
        self.dirlist = dirlist
        
        
#        self.list_amps = [i.split('Mesh',1)[1] for i in self.dirlist]
#        self.list_amps = [i.split('_',1)[0] for i in self.list_amps]
#        
#        
#        self.list_widths = [i.split('_',2)[1] for i in self.dirlist]
#        
#        self.df_amps_widths = pd.DataFrame(list(zip(self.list_amps,self.list_widths)),columns=['Amplitudes','Widths']) 
#        self.df_amps_widths =  self.df_amps_widths.sort_values(by=['Amplitudes','Widths'])
         
      
        # create fodername of the run for 2d:
        if self.dimensions == str('2'):
                self.run_folder = 'Mesh{:}_{:}_{:}_{:}'.format(
                       self.bump_amplitude,
                       self.bump_distribution_x,
                       self.bump_distribution_y,
                       self.prop) 
                
         # create fodername of the run for 2d:
        elif self.dimensions == None: 
                self.run_folder = 'Mesh{:}_{:}{:}_{:}'.format(
                       self.bump_amplitude,
                       self.bump_distribution_x,
                       self.bump_distribution_y,
                       self.prop)         
        
        
        
        # path to the directory of the model run:
        if self.dimensions is None:
            self.run_directory = os.path.join(self.res_folder,
                                           self.run_folder,
                                           self.sub_mesh_directory) 
    
        elif self.dimensions==str("2"):
             self.run_directory = os.path.join(self.res_folder,
                                           self.run_folder + '/') 
         
    
       
        # sort timesteps after execution (2d or 3d valid)
        if self.dimensions is None:
                self.dic_timesteps = glob.glob(self.run_directory + '*pvtu')
                
                
        elif self.dimensions==str("2"):
                self.dic_timesteps = glob.glob(self.run_directory + '*pvtu')
                
        self.dic_timesteps.sort(key=os.path.getmtime)
            
        # Get number of timesteps
        # while True:
        #    try:
        self.f_name = self.dic_timesteps[self.timestep]
            
        #    except IndexError:
        #        print('Given value for timestep is not valid. There is probably '
        #                  'no folder with given input parameters. Check directory.')
                
        self.num_timesteps = len(self.dic_timesteps)  
        
        print("Timesteps = ", self.num_timesteps)
        print(self.f_name)
           
        self.path_timestep = os.path.join(self.run_directory, self.f_name)
        self.dict_var = {}
        
        
        # Define xmlreader (file-input for VTU-path)
        # Create vtu objects, to collect properties, strings, variables etc
        # choose pvtu or vtu (2d or 3d)
  
        self.xmlReader = vtk.vtkXMLPUnstructuredGridReader()                            
        self.xmlReader.SetFileName(self.f_name)
        self.xmlReader.Update()
        
        # Get outputs of variables and geometry (2D and 3D)
        self.npoints = self.xmlReader.GetOutput().GetNumberOfPoints()
        self.ncells = self.xmlReader.GetOutput().GetNumberOfCells()
        self.bounds = self.xmlReader.GetOutput().GetBounds()
        self.narrays = self.xmlReader.GetOutput().GetPointData().GetNumberOfArrays()
        self.Points = vtk_to_numpy(self.xmlReader.GetOutput().GetPoints().GetData())
        self.Cells =  vtk_to_numpy(self.xmlReader.GetOutput().GetCells().GetData())
        
        self.x = self.Points[:,0]
        self.y = self.Points[:,1]
        self.z = self.Points[:,2]
        
        
        
        # Get outputs specified for 2D - case
        if self.dimensions==str("2"):
             self.sxy = vtk_to_numpy(self.xmlReader.GetOutput().GetPointData().GetArray(' sxy'))
             self.fs_upper = vtk_to_numpy(self.xmlReader.GetOutput().GetPointData().GetArray('fs upper'))
             self.fs_lower = vtk_to_numpy(self.xmlReader.GetOutput().GetPointData().GetArray('fs lower'))
             self.get_original_width = int(round(np.sqrt(self.bump_distribution_x)*2)*2)
        
        # Check for variables in list
        self.list_var_names = []
        for i in range(self.narrays):            
             self.name = self.xmlReader.GetPointArrayName(i)
             self.list_var_names.append(self.name)
            
            
            # Dictionary for variables
             self.dict_var_names = {}
             self.dict_var_values = {}
            
             for i in range(self.narrays):
        
                self.dict_var_names[i] = self.xmlReader.GetPointArrayName(i)
                self.dict_var_values[i] = self.xmlReader.GetOutput().GetPointData().GetArray(i)
                self.dict_var_names[self.xmlReader.GetPointArrayName(i)] \
                = self.xmlReader.GetOutput().GetPointData().GetArray(i)
        
        
        
    def get_scalar(self, 
                   scalar_name,
                   ):
         
         '''
          This method return the array or scalar as a numpy array. Check
          ModelRun().list_var_names for possible input.
          
         '''   
         
         self.scalar_name = scalar_name
         
         self.scalar = vtk_to_numpy(self.xmlReader.GetOutput().GetPointData().GetArray(str(self.scalar_name)))
         

        
         return self.scalar
     
     
        
      
    def get_scalar_real(self,
                        scalar_name
                        ):
        '''
        
        '''
        
        self.scalar_name = scalar_name
         
        self.scalar = vtk_to_numpy(self.xmlReader.GetOutput().GetPointData().GetArray(str(self.scalar_name)))
        
        
        neg = 2687
        pos = 2688       
        y_cells = 16
        x_half_cells = 167
        
        # Split domain to filter ghost cells
        x_neg = self.x[:neg]
        x_pos = self.x[pos:]
        y_neg = self.y[:neg]
        y_pos = self.y[pos:]
        
        x_neg_ind = np.argwhere(x_neg>=0)
        x_neg_new = np.delete(x_neg,x_neg_ind)
        y_neg_new = np.delete(y_neg,x_neg_ind)
        
        
        x_neg_mat = x_neg_new.reshape(y_cells,x_half_cells)
        x_pos_mat = x_pos.reshape(y_cells,x_half_cells)
        
        y_neg_mat = y_neg_new.reshape(y_cells,x_half_cells)
        y_pos_mat = y_pos.reshape(y_cells,x_half_cells)
         

        # Return
        x_mat = np.append(x_neg_mat,x_pos_mat,axis=1)
        y_mat = np.append(y_neg_mat,y_pos_mat,axis=1)
        
        if self.scalar.shape == self.x.shape:
            
            sxy_neg = self.scalar[:neg]
            sxy_pos = self.scalar[pos:]
            
            # Deleted for neg            
            sxy_neg_new = np.delete(sxy_neg,x_neg_ind)
            
            sxy_neg_mat = sxy_neg_new.reshape(y_cells,x_half_cells)
            sxy_pos_mat = sxy_pos.reshape(y_cells,x_half_cells)
            
            
            # Return
            
            scalar_mat = np.append(sxy_neg_mat,sxy_pos_mat,axis=1)
            
            return x_mat,y_mat,scalar_mat
        
        elif self.scalar.shape == self.Points.shape:
            
            self.scalar_x = self.scalar[:,0]
            self.scalar_y = self.scalar[:,1]
            self.scalar_z = self.scalar[:,2]
            
            sxy_neg_x = self.scalar_x[:neg]
            sxy_pos_x = self.scalar_x[pos:]
            
            sxy_neg_y = self.scalar_y[:neg]
            sxy_pos_y = self.scalar_y[pos:]
            
            sxy_neg_z = self.scalar_z[:neg]
            sxy_pos_z = self.scalar_z[pos:]
            
            
            # Deleted for neg            
            sxy_neg_new_x = np.delete(sxy_neg_x,x_neg_ind)
            sxy_neg_new_y = np.delete(sxy_neg_y,x_neg_ind)
            sxy_neg_new_z = np.delete(sxy_neg_z,x_neg_ind)
            
            sxy_neg_mat_x = sxy_neg_new_x.reshape(y_cells,x_half_cells)
            sxy_pos_mat_x = sxy_pos_x.reshape(y_cells,x_half_cells)
            
            sxy_neg_mat_y = sxy_neg_new_y.reshape(y_cells,x_half_cells)
            sxy_pos_mat_y = sxy_pos_y.reshape(y_cells,x_half_cells)
            
            sxy_neg_mat_z = sxy_neg_new_z.reshape(y_cells,x_half_cells)
            sxy_pos_mat_z = sxy_pos_z.reshape(y_cells,x_half_cells)
            
            # Return            
            scalar_mat_x = np.append(sxy_neg_mat_x,sxy_pos_mat_x,axis=1)
            scalar_mat_y = np.append(sxy_neg_mat_y,sxy_pos_mat_y,axis=1)
            scalar_mat_z = np.append(sxy_neg_mat_z,sxy_pos_mat_z,axis=1)
            
            return x_mat,y_mat,scalar_mat_x,scalar_mat_y,scalar_mat_z
            
            
            
            
        
    def cutter(self ,
               GL=None):
        
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
            self.GL = 1056000
        elif GL is not None:
            self.GL = GL
        
        
        
        # create a plane to cut; here it cuts in the XZ direction
        #(xz normal=(1,0,0);XY =(0,0,1),YZ =(0,1,0)
        self.plane_shelf =vtk.vtkPlane()
        self.plane_shelf.SetOrigin(self.GL,0,0)
        self.plane_shelf.SetNormal(1,0,0)

        #create cutter
        self.clipData_shelf = vtk.vtkClipDataSet()   
        self.clipData_shelf.SetClipFunction(self.plane_shelf)
        self.clipData_shelf.SetInputConnection(self.xmlReader.GetOutputPort())
        self.clipData_shelf.Update()
        
        # create plane to cut; end of shelf is cut due to boundary and mesh-
        # refinement adjustment for plotting
        self.plane_shelf_end = vtk.vtkPlane()
        self.plane_shelf_end.SetOrigin(1079000,0,0)
        self.plane_shelf_end.SetNormal(-1,0,0)
        
        #create cutter
        self.clipData = vtk.vtkClipDataSet()   
        self.clipData.SetClipFunction(self.plane_shelf_end)
        self.clipData.SetInputConnection(self.clipData_shelf.GetOutputPort())
        self.clipData.Update()
            
        
        return self.clipData
        

    
        
    def compute_convexhull(self):
        
        """
        Method: convexhull 
        ----------
        Computes boundary hull of cutted domain. Only valid if boundaries
        are rectangular (e.g. hull of 3D- Model)
        
        Parameters
        ----------
        -
        
        Returns:
        ----------
        Tuple(float): Points of convex-hull (x, y, z)
                      
        """
        
        
        self.triangulation = vtk.vtkDelaunay3D()
        self.triangulation.SetInputData(self.cutter().GetOutput())
        self.triangulation.Update()
        
        self.convexhull = vtk.vtkDataSetSurfaceFilter() 
        self.convexhull.SetInputConnection(self.triangulation.GetOutputPort())
        self.convexhull.Update()
        
        self.ch_points = self.convexhull.GetOutput()
        
        return self.ch_points
              

 
                  
    def selectMinz(self, x, y, z):
        
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
    

    
    def selectMaxz(self,x ,y, z):
        
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
        
        #======================================================================        
        # 3D  
        #====================================================================== 
        
        if self.dimensions is None:
            
            # Paraneter setup for calculation of hydrostatic thickness
            # !!!!!!!! Must be changed if input.sif file is changed!!!!!!!!
            self.p_w = 1000.0 # kg m−3 )
            self.p_i = 900.0  # ice 
            self.p_a = 2.0    # air (ρa =2kgm−3)
            self.H_a = 0      # firn corrections      
            
            #self.cutter = self.clipData
            self.clipped_area = self.cutter()
            self.zs = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray('zs'))
            self.zb = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray('zb'))
            self.points = vtk_to_numpy(self.clipped_area.GetOutput().GetPoints().GetData())
            
            
            # Get Points
            self.x = self.points[:,0]
            self.y = self.points[:,1]
            self.z = self.points[:,2]
    
            # Function to find min-values of z_b       
            self.zmin = self.selectMinz(self.x,self.y,self.z)
            self.zmax = self.selectMaxz(self.x,self.y,self.z)
    
            zs = self.zmax[1]
            zs = np.asarray(zs)
            xy_zmax = self.zmax[0]
           
            
            zb = self.zmin[1]
            zb = np.asarray(zb)
            xy_zmin = self.zmin[0]
            
    
            self.zb_new = zb
            self.zs_new = zs
            
            self.x_corr = xy_zmax[0] 
            self.y_corr = xy_zmax[1]
            
            # Calculation of real hydrostatic thickniss, modelled thickness and
            # deviation of hydrostatic equilibrium
            self.thick_model = -self.zs_new+self.zb_new       
            self.thick_calc = np.divide((self.p_w*self.zs_new),(self.p_w-self.p_i))                
            self.h_thickness = self.thick_model + self.thick_calc 
            
            
            
            # Return values of coordinates and calculated hydrostatic thickness
            return self.x_corr, self.y_corr, self.h_thickness, self.thick_model
        
        
        #======================================================================        
        # 2D 
        #====================================================================== 
        
        elif self.dimensions is not None:
            
            # Paraneter setup for calculation of hydrostatic thickness
            # !!!!!!!!!! Must be changed if input.sif file is changed!!!!!!!!
            
            self.p_w_2d = 1025.0 # kg m−3 )
            self.p_i_2d = 910.0  # ice (ρi = 918kgm−3)
            
            # pick surface and bottom coordinates
            self.fs_upper = vtk_to_numpy(self.xmlReader.GetOutput(). \
                                         GetPointData().GetArray('fs upper'))
            self.fs_lower = vtk_to_numpy(self.xmlReader.GetOutput(). \
                                         GetPointData().GetArray('fs lower'))
            self.points = vtk_to_numpy(self.xmlReader.GetOutput(). \
                                       GetPoints().GetData())
            
            
            self.x = self.points[:,0]
            self.y = self.points[:,1]
            
            # Rearrange matrix to fit input shape ()
            self.points = np.delete(self.points, 2, 1)
           
            # Group
            
            ind_fs_lower = np.where(self.fs_lower<0)
            ind_fs_upper = np.where(self.fs_upper>0)
            
            corr_fs_lower = self.fs_lower[ind_fs_lower] 
            corr_fs_upper = self.fs_upper[ind_fs_upper]
            x_new_upper =  self.x[ind_fs_upper]
            x_new_lower =  self.x[ind_fs_lower] 
            
            self.upper_model = self.y[ind_fs_upper]
            self.lower_model = self.y[ind_fs_lower] 
            
            # Real model thickness
            self.thick_model_2d = -corr_fs_upper+corr_fs_lower
            
            # Calculated thickness and thickness below sea-level
            self.thick_calc_2d = np.divide((self.p_w_2d*corr_fs_upper),(self.p_w_2d-self.p_i_2d)) 
            self.thick_calc_2d = -1 * self.thick_calc_2d

            # Final calculated hydrostatic thickness                    
            self.thick_calc_2d_bs = self.thick_calc_2d  + self.upper_model
            
            
            self.h_thickness = self.thick_model_2d + self.thick_calc_2d
            
            
            # Calculate peak deviation (maximum of modelled and calculated ht)
            self.deviation_list = []
            
            ht_calc = self.thick_calc_2d_bs
            ht_model = self.lower_model
            
            ht_calc_ind = int(abs(ht_calc.size/2)) 
            ht_model_ind = int(abs(ht_model.size/2))
        
            peak_calc = abs(ht_calc[ht_calc_ind])
            peak_model = abs(ht_model[ht_model_ind])
            
            self.deviation_list = 100-(100/peak_calc*peak_model)
            
            
            # Channel abutsment
            hydrostatic_deviation = self.thick_calc_2d_bs - self.lower_model
            
            peak = np.where(self.thick_calc_2d_bs==self.thick_calc_2d_bs.min())
            hd_peak = hydrostatic_deviation[peak]
            self.hd_abut = 100-(100/self.thick_calc_2d_bs[peak]*self.lower_model[peak])
            
            # Get indices of flank 
            hd_flank = np.max(hydrostatic_deviation)
            hd_flank_ind = np.where(hydrostatic_deviation == np.max(hydrostatic_deviation))
            
            self.deviation_list_flank = abs(100-(100/self.thick_calc_2d_bs[hd_flank_ind]*self.lower_model[hd_flank_ind]))
            
            return x_new_lower, self.thick_calc_2d_bs,self.upper_model, \
                self.lower_model, self.points, self.fs_lower, self.fs_upper, \
                self.deviation_list,self.deviation_list,self.hd_abut, -peak_model
     


    #@compute_hydrostatic_thickness               
    def compute_max_peak_deviation(self):
        pass
    
        '''
        Method not used yet (April,2019)
        
        '''
        
        # Calculate peak deviation (maximum of modelled and calculated ht)
        ht_calc = self.thick_calc_2d_bs
        ht_model = self.lower_model
        
        ht_calc_ind = int(abs(ht_calc.size/2)) 
        ht_model_ind = int(abs(ht_model.size/2))
        
        peak_calc = abs(ht_calc[ht_calc_ind])
        peak_model = abs(ht_model[ht_model_ind])
        
        # absolute difference (peak difference in m)
        abs_diff = peak_calc-peak_model
        
        return abs_diff
    
        
        
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
        self.levels = 19
        
        
        # Arange plane to cut
        self.line=vtk.vtkPlane()
        self.line.SetOrigin(self.xcoord,0,0)
        self.line.SetNormal(1,0,0)
        
        self.clipped_area_3D = self.cutter()
        
        # Get cutter
        self.cutter_3D = vtk.vtkCutter()   
        self.cutter_3D.SetCutFunction(self.line)
        self.cutter_3D.SetInputConnection(self.clipped_area_3D.GetOutputPort())
        self.cutter_3D.Update()
        self.line_points = vtk_to_numpy(self.cutter_3D.GetOutput().GetPoints().GetData())
        
        # Write numpy array
        self.out = vtk_to_numpy(self.cutter_3D.GetOutput().GetPointData().GetArray(self.var))
        
        df = pd.DataFrame({'y':self.line_points[:,1],'z':self.line_points[:,2],'scalar':self.out})
        df = df.sort_values(by=['y','z'])
        
        
        mat_sorted = df.values
        
        y_sort = mat_sorted[:,0]
        z_sort = mat_sorted[:,1]
        sxy_sort = mat_sorted[:,2]
        
        self.y_mat = np.flipud(y_sort.reshape(int(y_sort.size/self.levels),self.levels).T)
        self.z_mat = np.flipud(z_sort.reshape(int(y_sort.size/self.levels),self.levels).T)
        self.scalar_mat = np.flipud(sxy_sort.reshape(int(y_sort.size/self.levels),self.levels).T)        
        


        return self.out, self.line_points, self.y_mat, self.z_mat, self.scalar_mat


        
 
    def compute_concavehull(self,
                            xcoord):
        
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
        
        #Clipped area
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
        
        # Upper and lower boundary
        lower_boundary = grouped.min()
        upper_boundary = grouped.max()
        
        # Lower boundary for calculated hydrostatic thickness
        original_thickness = grouped_calc.min()
        
        
        self.hull_panda = upper_boundary,lower_boundary,original_thickness
        
        
            
               
        
        return upper_boundary, lower_boundary, original_thickness
    
    
    
    
    def compute_melt_rate(self):
        pass
        
        self.a_smb = 0.3 # m/a
        self.clipped_area = self.cutter()
        self.zs = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray('zs'))
        self.u  = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray('velocity'))
        self.ux = self.u[:,0]
        print(self.ux.shape)
        print(self.zs.shape)
        self.melt = self.a_smb-((self.zs*self.ux)/(1-(self.p_i/self.p_w)))
        
        
        return self.melt
    
    
    def compute_second_invariant(self):
        
        '''
        '''
        
        
        if self.dimensions is None:
            
            self.clipped_area = self.cutter()
            self.sxx = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray('sxx'))
            self.syy = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray(' syy'))
            self.szz = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray('szz'))
        
            self.j_2 = 1/2*((self.sxx**2)+(self.syy**2)+(self.szz**2))
            
            
            # pick surface and bottom coordinates
            self.fs_upper = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray('zs'))
    
            self.points = vtk_to_numpy(self.clipped_area.GetOutput().GetPoints().GetData())
            self.ind_fs = np.where(self.fs_upper!=0)
            
            self.x = self.points[:,0]
            self.y = self.points[:,1]
            
            self.x_upper = self.x[self.ind_fs]
            self.y_upper = self.y[self.ind_fs]
            self.j_2_upper = self.j_2[self.ind_fs]
            
             
             
             
        
        return self.x_upper,self.y_upper,self.j_2_upper
    
    
    
    
    
    
    def get_top_scalar(self,scalar_name):
        
        '''
        '''
        
       
        
        self.clipped_area = self.cutter()
        self.scalar = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray(str(scalar_name)))
        
        if scalar_name == str('velocity'):
            self.scalar = self.scalar[:,1]
            
        
        # pick surface and bottom coordinates
        self.fs_upper = vtk_to_numpy(self.clipped_area.GetOutput().GetPointData().GetArray('zs'))

        self.points = vtk_to_numpy(self.clipped_area.GetOutput().GetPoints().GetData())
        self.ind_fs = np.where(self.fs_upper!=0)
        
        self.x = self.points[:,0]
        self.y = self.points[:,1]
        
        self.x_upper = self.x[self.ind_fs]
        self.y_upper = self.y[self.ind_fs]
        self.scalar_top = self.scalar[self.ind_fs]
        
             
             
             
        
        return self.x_upper,self.y_upper,self.scalar_top
        
        
        
    
    
  
        
        
         
      
        
        
         
         
     
   
           
           
           
                        
    
            
            
    

