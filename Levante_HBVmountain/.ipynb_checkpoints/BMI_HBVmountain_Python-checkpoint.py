#!/usr/bin/env python
# coding: utf-8

# In[14]:
import julia
from julia import Main
from bmipy import Bmi
from pandas import read_csv
import numpy as np
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import datetime


# In[22]:
Main.include("Refactoring/structs.jl")
Main.include('Refactoring/Julia_BMI_model.jl')
Main.include("Refactoring/parameterselection.jl")
Main.include("Refactoring/processes_buckets.jl")
Main.include("Refactoring/elevations.jl")
Main.include("Refactoring/allHRU.jl")
Main.include("Refactoring/run_model.jl");
Main.include("Refactoring/Config_model.jl")
Main.include("Refactoring/Finalize_Model.jl")
Main.include("Refactoring/Potential_Evaporation.jl")
Main.include("Refactoring/Preprocessing.jl")
Main.include("Refactoring/Units.jl")

build_HBVmountain_model = Main.build_HBVmountain_model
setup = Main.setup
initialize = Main.initialize
update = Main.update
update_until = Main.update_until
get_component_name = Main.get_component_name
get_input_var_names = Main.get_input_var_names
get_output_var_names = Main.get_output_var_names
get_start_time = Main.get_start_time
get_current_time = Main.get_current_time
get_end_time = Main.get_end_time
get_time_step = Main.get_time_step
get_time_units = Main.get_time_units
get_var_type = Main.get_var_type
get_var_units = Main.get_var_units
get_var_itemsize = Main.get_var_itemsize
get_var_nbytes = Main.get_var_nbytes
get_value = Main.get_value
get_value_at_indices = Main.get_value_at_indices
get_value_ptr = Main.get_value_ptr
get_var_location = Main.get_var_location
get_input_item_count = Main.get_input_item_count
get_output_item_count = Main.get_output_item_count
set_value = Main.set_value
set_value_at_indices = Main.set_value_at_indices

Parameters = Main.Parameters
Slow_Paramters = Main.Slow_Paramters
# In[22]:
class BMI_HBVmountain(Bmi):
    def __init__(self):
        self.model = build_HBVmountain_model() 
        
    def setup(self):
        return setup()
 
    def initialize(self, config_file):
        initialize(self.model, config_file);

    def update(self):
        update(self.model)
    
    def update_until(self, time):
        update_until(self.model, time)
        
    def finalize(self):
        self.model = None
    # info
    def get_component_name(self):
        return get_component_name(self.model)

    def get_input_var_names(self):
        return Main.get_input_var_names(self.model)

    def get_output_var_names(self):
        return get_output_var_names(self.model)

    # time
    def get_start_time(self):
        return get_start_time(self.model)

    def get_current_time(self):
        return get_current_time(self.model)

    def get_end_time(self):
        return get_end_time(self.model)

    def get_time_step(self):
        return get_time_step(self.model)

    def get_time_units(self):
        return get_time_units(self.model)

    # vars
    def get_var_type(self, name):
        return get_var_type(self.model, name)

    def get_var_units(self, name):
        return get_var_units(self.model, name)

    def get_var_itemsize(self, name):
        return get_var_itemsize(self.model, name)

    def get_var_nbytes(self, name):
        return get_var_nbytes(self.model, name)

    # getter
    def get_value(self, name, dest):
        return get_value(self.model, name, dest)


    def get_value_at_indices(self, name, dest, inds):
        return get_value_at_indices(self.model, name, dest, inds)


    def get_value_ptr(self, name):
        return get_value_ptr(self.model, name)

    def get_var_location(self, name):
        return get_var_location(self.model, name)
    
    def get_input_item_count(self):
        return get_input_item_count(self.model)
    
    def get_output_item_count(self):
        return get_output_item_count(self.model)    
    
    # setter
    def set_value(self, name, value):
        set_value(self.model, name, value)


    def set_value_at_indices(self, name, inds, value):
        set_value_at_indices(self.model, name, inds, value)
    
    # grid
    def get_var_grid(self, name):
        return Main.get_var_grid(self.model, name)
    
    def get_grid_rank(self, grid):
        return Main.get_grid_rank(self.model, grid)

    def get_grid_size(self, grid):
        return Main.get_grid_size(self.model, grid)

    def get_grid_type(self, grid):
        return Main.get_grid_type(self.model, grid)

    def get_grid_shape(self, grid):
        return Main.get_grid_shape(self.model, grid)

    def get_grid_x(self, grid, x):
        return Main.get_grid_x(self.model, grid, x)

    def get_grid_y(self, grid, y):
        return Main.get_grid_x(self.model, grid, y)

    def get_grid_z(self, grid, z):
        return Main.get_grid_x(self.model, grid, z)

    def get_grid_spacing(self, grid, spacing):
        return Main.get_grid_spacing(self.model, grid, spacing)

    def get_grid_origin(self, grid, origin):
        return Main.get_grid_origin(self.model, grid, origin)

    def get_grid_node_count(self, grid):
        return Main.get_grid_node_count(self.model, grid)

    def get_grid_edge_count(self, grid):
        return Main.get_grid_edge_count(self.model, grid)

    def get_grid_face_count(self, grid):
        return Main.get_grid_face_count(self.model, grid)

    def get_grid_edge_nodes(self, grid, edge_nodes):
        return Main.get_grid_edge_nodes(self.model, grid, edge_nodes)
    
    def get_grid_face_edges(self, grid, face_edges):
        return Main.get_grid_face_edges(self.model, grid, face_edges)
       
    def get_grid_face_nodes(self, grid, face_nodes):
        return Main.get_grid_face_nodes(self.model, grid, face_nodes)

    def get_grid_nodes_per_face(self, grid, nodes_per_face):
        return Main.get_grid_nodes_per_face(self.model, grid, nodes_per_face)





