# -*- coding: utf-8 -*-
"""
This module holds the Timestep class
"""

class timestep():
    """
    Timestep class for a vtu file of one timestep in a model run
    
    Attributes
    ----------
        data
    """
    def __init__(self, filename):
        """
        contructor
        """
        self.filepath = filename
        self.data = {}
        
    def load_vtu(self):
        """
        load vtu file
        """
        # einlade code hier
        # sollte am ende data auswerfen
        data = {'variable1': data1,
                'variable2': data2}
        return self.data
    
    def plot(self, variable):
        
        plot_data = copy(self.data[variable])
          # def timestep(self,timestep):
            
