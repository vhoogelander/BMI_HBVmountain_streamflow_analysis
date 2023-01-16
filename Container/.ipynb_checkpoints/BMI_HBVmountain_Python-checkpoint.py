#!/usr/bin/env python
# coding: utf-8
# In[14]:
from utils import *


class BMI_HBVmountain(Bmi):
    """
    Creates BMI HBV-mountain model object
    """
    def __init__(self, forcing_netcdf=None, path_to_shapefile=None, path_to_dem=None, path_to_nlcd=None):

        self.model = build_HBVmountain_model() #Julia function: create model object
        self.forcing_netcdf = forcing_netcdf
        self.path_to_shapefile = path_to_shapefile
        self.path_to_dem = path_to_dem
        self.path_to_nlcd = path_to_nlcd

    def setup(self, forcing_netcdf=None, path_to_shapefile=None, path_to_dem=None, path_to_nlcd=None, Discharge=0.0, Total_Evaporation=0.0, Snow_Extend=np.zeros(1),
              bare_storage=Storages(0, np.zeros(4), np.zeros(4), np.zeros(4), 0), forest_storage=Storages(0, np.zeros(4), np.zeros(4), np.zeros(4), 0), 
              grass_storage=Storages(0, np.zeros(4), np.zeros(4), np.zeros(4), 0), rip_storage=Storages(0, np.zeros(4), np.zeros(4), np.zeros(4), 0), Slowstorage=0.0, 
              Waterbalance=0.0, Glacier=[0,0,0,0], Area=0,
              bare_parameters=None, forest_parameters=None, grass_parameters=None, rip_parameters=None, slow_parameters=None, 
              Temp_Thresh=0.0, Meltfactor=3.0, Mm=0.01, Ratio_Pref=0.4, Kf=1.0, Ks=0.05, Ce=0.6, Soilstoragecapacity_Bare=40.0, beta_Bare=1.0,
              Interceptioncapacity_Forest=7.0, Soilstoragecapacity_Forest=500.0, beta_Forest=1.0, Interceptioncapacity_Grass=3.0, Soilstoragecapacity_Grass=250.0, 
              beta_Grass=1.0, Interceptioncapacity_Rip=4.0, Soilstoragecapacity_Rip=100.0, beta_Rip=1.0, Kf_Rip=1.3, Ratio_Riparian=0.1, 
              Elevation=Elevations(100,100,500,250,250),Total_Elevationbands=4, Precipitation_gradient=0, Elevation_Percentage=[0.25,0.25,0.25,0.25], 
              bare_input=HRU_Input([0,0,0,0], 0, np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.0, np.zeros(4), 0, 0.0), forest_input=HRU_Input([0.3,0.3,0.3,0.1],
              0.5,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0), 
              grass_input=HRU_Input([0.1,0.3,0.3,0.3], 0.47,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0), 
              rip_input=HRU_Input([0.0,0.0,0.3,0.7], 0.03,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0), 
              Precipitation=np.zeros((0,0)), Temperature=np.zeros((0,0)), ETP=np.zeros(0), Date=None, Current_Date=None, 
              Sunhours = [8.87, 10.30, 11.88, 13.65, 15.13, 15.97, 15.58, 14.25, 12.62, 11.87, 9.28, 8.45],  Units=HBVmountain_units()):
        """
        Returns a configuration file for the setup of BMI HBV-mountain model. Preprocessing of model settings is done using the catchment's shapefile, DEM raster file and NLCD landuse (Must be in WGS84)
        Parameters
        ----------
        forcing_netcdf: netCDF4.Dataset
        path_to_shapefile: String path catchment shapefile. Must be downloaded from GRDC (WGS84)
        path_to_dem: String path to rasterfile (WGS84)
        path_to_nlcd: String path to rasterfile (WGS84)
        Discharge : Float
        Total_Evaporation : Float
        Snow_Extend : Array with length Total_Elevationbands
        bare_storage : HBV-mountain storage object
        forest_storage : HBV-mountain storage object
        grass_storage : HBV-mountain storage object
        rip_storage : HBV-mountain storage object
        Slowstorage : Float
        Waterbalance : Float
        Glacier : Array with length Total_Elevationbands
        Area : 
        bare_parameters : HBV-mountain parameter object
        forest_parameters : HBV-mountain parameter object
        grass_parameters : HBV-mountain parameter object
        rip_parameters : HBV-mountain parameter object
        slow_parameters : HBV-mountain slow parameter object
        Elevation : HBV-mountain elevations object
        Total_Elevationbands : Int
        Precipitation_gradient : Float
        Elevation_Percentage : 
        bare_input : HBV-mountain HRU_Input object
        forest_input : HBV-mountain HRU_Input object
        grass_input : HBV-mountain HRU_Input object
        rip_input : HBV-mountain HRU_Input object
        Precipitation : Array (len data x Total_Elevationbands) 
        Temperature : Array (len data x Total_Elevationbands) 
        ETP : Array 
        Date : Datetime
        Current_Date : Datetime
        Sunhours : Array
        Units : HBVmountain_model_units object            
        """
        
        
        if bare_parameters == None:
            bare_parameters = Parameters(beta_Bare, Ce, 0, 0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoragecapacity_Bare, Temp_Thresh)
        if forest_parameters == None:
            forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoragecapacity_Forest, Temp_Thresh)
        if grass_parameters == None:
            grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoragecapacity_Grass, Temp_Thresh)
        if rip_parameters == None:
            rip_parameters = Parameters(beta_Rip, Ce, 0, Interceptioncapacity_Rip, Kf_Rip, Meltfactor, Mm, Ratio_Pref, Soilstoragecapacity_Rip, Temp_Thresh)
        if slow_parameters == None:
            slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)
        
        if self.forcing_netcdf != None:
            forcing = generate_forcing_from_NETCDF(self.forcing_netcdf)

            Date = forcing.index.values
            Current_Date = forcing.index.values[0]
            Precipitation = forcing.prec.values.reshape(len(forcing),1)
            Temperature = forcing.temp.values.reshape(len(forcing),1)
        
        
        if self.path_to_shapefile != None:
            shapefile = gpd.read_file(self.path_to_shapefile)
            lon, lat = shapefile.iloc[0].geometry.centroid.x, shapefile.iloc[0].geometry.centroid.y
            Sunhours = get_monthly_sunhours(lat, lon)
            
            landuse_elevation = generate_landuse_per_elevation(str_path_to_elevation=self.path_to_dem, str_path_to_landuse=self.path_to_nlcd, str_path_to_shapefile=self.path_to_shapefile, nbands=4)
            elevation_step = landuse_elevation[0][2]
            min_elevation, max_elevation, mean_elevation = landuse_elevation[0][3], landuse_elevation[0][4], landuse_elevation[0][5]
            Elevation = Elevations(elevation_step,min_elevation,max_elevation,mean_elevation,mean_elevation)
            Glacier = landuse_elevation[1][0]
            
            bare_input = HRU_Input(landuse_elevation[2], landuse_elevation[0][1][0], landuse_elevation[1][1], [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0)
            forest_input = HRU_Input(landuse_elevation[3], landuse_elevation[0][1][1], np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0)
            grass_input = HRU_Input(landuse_elevation[4], landuse_elevation[0][1][2], np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0)
            rip_input = HRU_Input(landuse_elevation[5], landuse_elevation[0][1][3], np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0)
            Total_Elevationbands = len(landuse_elevation[0][0])
            Elevation_Percentage = landuse_elevation[0][0]
        
        return setup(Discharge, Total_Evaporation, Snow_Extend, bare_storage,forest_storage, grass_storage, rip_storage,
        Slowstorage, Waterbalance, Glacier, Area, bare_parameters, forest_parameters, grass_parameters, rip_parameters, slow_parameters,
        Elevation, Total_Elevationbands, Precipitation_gradient, Elevation_Percentage, bare_input, forest_input, grass_input, rip_input, Precipitation, Temperature, ETP, Date, Current_Date, Sunhours, Units)
 
    def initialize(self, config_file):
        """
        Initializes BMI HBV-mountain model using a configuration file
        """
        initialize(self.model, config_file);

    def update(self):
        """
        Updates the BMI HBV-mountain model for one timestep
        """
        update(self.model)
    
    def update_until(self, time):
        """
        Updates the BMI HBV-mountain model for until given time. Time must be a string
        """
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





