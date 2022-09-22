using Dates

mutable struct HRU_Input
    #inputs (alphabetic order)
    Area_Elevations::Array{Float64,1}
    Area_HRU:: Float64
    Area_Glacier::Array{Float64,1} # smaller than 1
    Elevation_Count::Array{Int64}
    Nr_Elevationbands:: Int64
    Catchment_Elevation::Tuple
    Snow_Redistribution::Tuple
    #Potential_Evaporation::Array{Float64,1} 
    Potential_Evaporation_Mean:: Float64
    Precipitation::Array{Float64,1}
    Riparian_Discharge:: Float64 #only necessary for riparian HRU
    Temp_Elevation::Array{Float64,1}
    Total_Effective_Precipitation::Float64
    Total_Interception_Evaporation::Float64
end

mutable struct Parameters
    # parameters (alphabetic order)
    beta:: Float64
    Ce:: Float64
    Drainagecapacity:: Float64 #only necessary for riparian HRU
    Interceptionstoragecapacity:: Float64
    Kf:: Float64
    Meltfactor:: Float64
    Mm:: Float64
    #Percolationcapacity:: Float64 #only necessary for hillslope HRU
    Ratio_Pref:: Float64 #only necessary for hillslope HRU
    Soilstoragecapacity:: Float64
    Temp_Thresh:: Float64
end

mutable struct Slow_Paramters
    Ks:: Float64
    Ratio_Riparian:: Float64
end

mutable struct Storages
    Fast:: Float64
    Interception::Array{Float64,1}
    Snow::Array{Float64,1}
    Snow_Cover::Array{Float64,1}
    Soil:: Float64
end

mutable struct Outflows
    Fast_Discharge:: Float64
    GWflow:: Float64 #only necessary for hillslope HRU
    Soil_Evaporation:: Float64
    Interception_Evaporation:: Float64
end

mutable struct Elevations
    Thickness_Band:: Float64
    Min_elevation:: Float64
    Max_elevation:: Float64
    Measured_Prec_Elevation:: Float64
    Measured_Temp_Elevation:: Float64
end

mutable struct Drought
    #inputs (alphabetic order)
    Nr_Drought_Days_Past::Array{Float64,1}
    Nr_Drought_Days_Future::Array{Float64,1}
    Nr_Drought_Events_Past::Array{Float64,1}
    Nr_Drought_Events_Future::Array{Float64,1}
    Max_Drought_Length_Past::Array{Float64,1}
    Max_Drought_Length_Future::Array{Float64,1}
    Mean_Drought_Length_Past::Array{Float64,1}
    Mean_Drought_Length_Future::Array{Float64,1}
    Max_Deficit_Past::Array{Float64,1}
    Max_Deficit_Future::Array{Float64,1}
    Mean_Deficit_Past::Array{Float64,1}
    Mean_Deficit_Future::Array{Float64,1}
    Total_Deficit_Past::Array{Float64,1}
    Total_Deficit_Future::Array{Float64,1}
    Max_Intensity_Past::Array{Float64,1}
    Max_Intensity_Future::Array{Float64,1}
    Mean_Intensity_Past::Array{Float64,1}
    Mean_Intensity_Future::Array{Float64,1}
end

mutable struct Drought_Extremes
    Longest_Drought_Length_Past::Array{Float64,1}
    Longest_Drought_Length_Future ::Array{Float64,1}
    Longest_Drought_Deficit_Past::Array{Float64,1}
    Longest_Drought_Deficit_Future::Array{Float64,1}
    Longest_Drought_Start_Past::Array{Float64,1}
    Longest_Drought_Start_Future::Array{Float64,1}
    Longest_Drought_End_Past::Array{Float64,1}
    Longest_Drought_End_Future::Array{Float64,1}
    Severest_Drought_Length_Past::Array{Float64,1}
    Severest_Drought_Length_Future ::Array{Float64,1}
    Severest_Drought_Deficit_Past::Array{Float64,1}
    Severest_Drought_Deficit_Future::Array{Float64,1}
    Severest_Drought_Start_Past::Array{Float64,1}
    Severest_Drought_Start_Future::Array{Float64,1}
    Severest_Drought_End_Past::Array{Float64,1}
    Severest_Drought_End_Future::Array{Float64,1}
end


mutable struct HBVmountain_model_units
    Discharge::Union{Nothing, String}
    Total_Evaporation::Union{Nothing, String}
    Snow_Extend::Union{Nothing, String}
    bare_storage::Union{Nothing, Array{String, 1}}
    forest_storage::Union{Nothing, Array{String, 1}}
    grass_storage::Union{Nothing, Array{String, 1}}
    rip_storage::Union{Nothing, Array{String, 1}}
    Slowstorage::Union{Nothing, String}
    Waterbalance::Union{Nothing, String}
    Glacier::Union{Nothing, String}
    #Model settings
    Area
    bare_parameters::Union{Nothing, Array{String, 1}}
    forest_parameters::Union{Nothing, Array{String, 1}}
    grass_parameters::Union{Nothing, Array{String, 1}}
    rip_parameters::Union{Nothing, Array{String, 1}}
    slow_parameters::Union{Nothing, Array{String, 1}}
    Elevation::Union{Nothing, Array{String, 1}}
    Total_Elevationbands ::Union{Nothing, String}
    Precipitation_gradient ::Union{Nothing, String}
    Elevation_Percentage 
    #Initial input settings 
    bare_input::Union{Nothing, Array{Union{Nothing, String}}}
    forest_input::Union{Nothing, Array{Union{Nothing, String}}}
    grass_input::Union{Nothing, Array{Union{Nothing, String}}}
    rip_input::Union{Nothing, Array{Union{Nothing, String}}}
    #Forcing
    Precipitation::Union{Nothing, String}
    Temperature::Union{Nothing, String}
    ETP::Union{Nothing, String}
    Date::Union{Nothing, String}
    Current_Date::Union{Nothing, String}
    Sunhours::Union{Nothing, String}
    
end

# module BasicModelInterface


mutable struct HBVmountain_model
    #Storages and discharge
    Discharge::Union{Nothing, Float64}
    Total_Evaporation::Union{Nothing, Float64}
    Snow_Extend::Union{Nothing, Array{Float64,1}}
    bare_storage::Union{Nothing, Storages}
    forest_storage::Union{Nothing, Storages}
    grass_storage::Union{Nothing, Storages}
    rip_storage::Union{Nothing, Storages}
    Slowstorage::Union{Nothing, Float64}
    Waterbalance::Union{Nothing, Float64}
    Glacier::Union{Nothing, Array{Float64,1}}
    #Model settings
    Area
    bare_parameters::Union{Nothing, Parameters}
    forest_parameters::Union{Nothing, Parameters}
    grass_parameters::Union{Nothing, Parameters}
    rip_parameters::Union{Nothing, Parameters}
    slow_parameters::Union{Nothing, Slow_Paramters}
    Elevation::Union{Nothing, Elevations}
    Total_Elevationbands ::Union{Nothing, Int64}
    Precipitation_gradient ::Union{Nothing, Float64}
    Elevation_Percentage 
    #Initial input settings
    bare_input::Union{Nothing, HRU_Input}
    forest_input::Union{Nothing, HRU_Input}
    grass_input::Union{Nothing, HRU_Input}
    rip_input::Union{Nothing, HRU_Input}
    #Forcing
    Precipitation::Union{Nothing, Array{Float64,2}}
    Temperature::Union{Nothing, Array{Float64,2}}
    ETP::Union{Nothing, Array{Float64,1}}
    #Other
    Date::Union{Nothing, Array{Date,1}}
    Current_Date::Union{Nothing, Date}
    Sunhours::Union{Nothing, Array{Float64, 1}}
    Units::Union{Nothing, HBVmountain_model_units}
end

