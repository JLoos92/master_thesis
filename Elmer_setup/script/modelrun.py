# -*- coding: utf-8 -*-

"""
this module holds the ModelRun class
"""
import os
from timestep import Timestep

class ModelRun:
    """
    This class holds the vtu files of the TimeSteps
    """
    
    def __init__(self, 
                 bump_amplitude, 
                 bump_distribution, 
                 bump_location, 
                 grid_refinement):
        """
        ModelRun constructor
        
        Parameters
        ----------
        bump_size : int
            length of the bump
        bump_amplitude : int
            height of the bump
        grid_refinement : bool
            whether the grid is refined or not
        """
        
        home_directory = '/Volumes/esd01/jloos/data_small'

        # change output directory according to mesh refinement:
        if grid_refinement:
            refinement_folder = os.path.join(home_directory,'refined')
        else:
            refinement_folder = os.path.join(home_directory,'uniform')
        # create fodername of the run:    
        run_folder = 'Mesh{:}_{:}_{:}'.format(bump_amplitude,bump_distribution,bump_location) 
        # path to the directory of the model run:
        self.run_directory = os.path.join(refinement_folder,run_folder)
        # construct Timestep objects for each timestep:
        self.timesteps = {}
        timesteps = os.listdir(run_directory)
        for step in range(len(timesteps)):
            self.timesteps[str(step)] = Timestep(os.path.join(self.run_directory,timesteps[step]))
            
    def load_timesteps(self, timesteps_to_consider='all'):
        """
        loads the vtu files indicated with steps
        """
        if timesteps_to_consider == 'all':
            self.timesteps.keys
        for step in timesteps_to_consider:
            self.timesteps[step].load_vtu()
    

    def plot_timesteps(self, timesteps_to_consider, variable):
        """
        plots the specified timesteps
        """
        for step in timesteps_to_consider:
            self.timesteps[step].plot(variable)
        
    

