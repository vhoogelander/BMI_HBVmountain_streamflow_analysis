function replace_HBVmountainmodel(model::HBVmountain_model, config_file::HBVmountain_model)
    model.Discharge = config_file.Discharge
    model.Snow_Extend = config_file.Snow_Extend
    model.bare_storage = config_file.bare_storage
    model.forest_storage = config_file.forest_storage
    model.grass_storage = config_file.grass_storage
    model.rip_storage = config_file.rip_storage
    model.Slowstorage = config_file.Slowstorage
    model.Waterbalance = config_file.Waterbalance
    model.Glacier = config_file.Glacier
    #Model settings
    model.Area = config_file.Area
    model.bare_parameters = config_file.bare_parameters
    model.forest_parameters = config_file.forest_parameters
    model.grass_parameters = config_file.grass_parameters
    model.rip_parameters = config_file.rip_parameters
    model.slow_parameters = config_file.slow_parameters
    model.Elevation = config_file.Elevation
    model.Total_Elevationbands = config_file.Total_Elevationbands
    model.Precipitation_gradient = config_file.Precipitation_gradient
    model.Elevation_Percentage =  config_file.Elevation_Percentage
    #Initial input settings
    model.bare_input = config_file.bare_input
    model.forest_input = config_file.forest_input
    model.grass_input = config_file.grass_input
    model.rip_input = config_file.rip_input
    #Forcing
    model.Precipitation = config_file.Precipitation
    model.Temperature = config_file.Temperature
    model.ETP = config_file.ETP
    model.Date = config_file.Date
    model.Current_Date = config_file.Current_Date
    model.Units = config_file.Units
end 
