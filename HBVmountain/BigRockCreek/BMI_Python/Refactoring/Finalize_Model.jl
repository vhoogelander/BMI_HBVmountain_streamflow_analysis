function finalize_HBVmountainmodel(model::HBVmountain_model)
    model.Discharge = nothing
    model.Snow_Extend = nothing
    model.bare_storage = nothing
    model.forest_storage = nothing
    model.grass_storage = nothing
    model.rip_storage = nothing
    model.Slowstorage = nothing
    model.Waterbalance = nothing
    model.Glacier = nothing
    #Model settings
    model.Area = nothing
    model.bare_parameters = nothing
    model.forest_parameters = nothing
    model.grass_parameters = nothing
    model.rip_parameters = nothing
    model.slow_parameters = nothing
    model.Elevation = nothing
    model.Total_Elevationbands = nothing
    model.Precipitation_gradient = nothing
    model.Elevation_Percentage =  nothing
    #Initial input settings
    model.bare_input = nothing
    model.forest_input = nothing
    model.grass_input = nothing
    model.rip_input = nothing
    #Forcing
    model.Precipitation = nothing
    model.Temperature = nothing
    model.ETP = nothing
    model.Date = nothing
    model.Current_Date = nothing
    model.Units = nothing
end 
