using DataFrames
using CSV
using Plots
using DocStringExtensions
using Dates
using NetCDF
using Statistics
using Pkg

function build_HBVmountain_model()
    return HBVmountain_model(nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    #Model settings
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    #Initial input settings
    nothing,
    nothing,
    nothing,
    nothing,
    #Forcing
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing)
end

# Create configuration file which can be used for initialization of the model
function setup(Discharge=0.0, Snow_Extend=zeros(1), bare_storage=Storages(0, zeros(1), zeros(1), zeros(1), 0), forest_storage=Storages(0, zeros(1), zeros(1), zeros(1), 0), 
        grass_storage=Storages(0, zeros(1), zeros(1), zeros(1), 0), rip_storage=Storages(0, zeros(1), zeros(1), zeros(1), 0), Slowstorage=0.0, Waterbalance=0.0, Glacier=zeros(1), Area=0, 
        bare_parameters=Parameters(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), forest_parameters=Parameters(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 
        grass_parameters=Parameters(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), rip_parameters=Parameters(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 
        slow_parameters=Slow_Paramters(0.0, 0.0), Elevation=Elevations(0,0,0,0,0),Total_Elevationbands=0, Precipitation_gradient=0, Elevation_Percentage=0, bare_input=HRU_Input([1], 0.0, [0], [1], 1, Tuple(0), Tuple(0), 0, [0], 0.0, [0], 0, 0.0), 
        forest_input=HRU_Input([1], 0.0, [0], [1], 1, Tuple(0), Tuple(0), 0, [0], 0.0, [0], 0, 0.0), 
        grass_input=HRU_Input([1], 0.0, [0], [1], 1, Tuple(0), Tuple(0), 0, [0], 0.0, [0], 0, 0.0), 
        rip_input=HRU_Input([1], 0.0, [0], [1], 1, Tuple(0), Tuple(0), 0, [0], 0.0, [0], 0, 0.0), 
        Precipitation=zeros(0,0), Temperature=zeros(0,0), ETP=zeros(0), Date=Date[], Current_Date=DateTime(0), Sunhours = [0], Units=HBVmountain_units())
    
    config_file = HBVmountain_model(Discharge, Snow_Extend, bare_storage,forest_storage, grass_storage, rip_storage,
        Slowstorage, Waterbalance, Glacier, Area, bare_parameters, forest_parameters, grass_parameters, rip_parameters, slow_parameters,
        Elevation, Total_Elevationbands, Precipitation_gradient, Elevation_Percentage, bare_input, forest_input, grass_input, rip_input, Precipitation, Temperature, ETP, Date, Current_Date, Sunhours, Units);
    return config_file
end

# uses function in Config_model.jl file
function initialize(model, config_file = nothing)
    if config_file != nothing
        replace_HBVmountainmodel(model, config_file) 
    else
        error("must give config_file"::AbstractString)
    end
end

# Updates model for one time step, uses the function in Refactoring/run_model 
function update(model::HBVmountain_model)
    #Update model for next timestep
    t = findall(x -> x == model.Current_Date, model.Date)[1]
    
    # Without Epot observations
    if sum(model.ETP) == 0
        model.ETP = getEpot_Daily_thornthwaite(model.Temperature[:,1], model.Date, model.Sunhours) 
    end
    
    #Without multiple Elevations
    if model.Elevation.Max_elevation == 0
        Discharge, Snow_Extend, Waterbalance, bare_storage, forest_storage, grass_storage, rip_storage, Slowstorage =  run_model_glacier(model.Area, model.ETP[t], model.Glacier, model.Precipitation[t,:], model.Temperature[t,:],
                        model.bare_input, model.forest_input, model.grass_input, model.rip_input,
                        model.bare_storage, model.forest_storage, model.grass_storage, model.rip_storage, model.Slowstorage,
                        model.bare_parameters, model.forest_parameters, model.grass_parameters, model.rip_parameters, model.slow_parameters, model.Total_Elevationbands, model.Elevation_Percentage)
    #With Elevations
    else
        Precipitation = getprecipitationatelevation(model.Elevation, model.Precipitation_gradient, model.Precipitation[t,:])[2][1,:]
        Temperature = gettemperatureatelevation(model.Elevation, model.Temperature[t,:])[2][1,:]
        Discharge, Snow_Extend, Waterbalance, bare_storage, forest_storage, grass_storage, rip_storage, Slowstorage =  run_model_glacier(model.Area, model.ETP[t], model.Glacier, Precipitation, Temperature,
                        model.bare_input, model.forest_input, model.grass_input, model.rip_input,
                        model.bare_storage, model.forest_storage, model.grass_storage, model.rip_storage, model.Slowstorage,
                        model.bare_parameters, model.forest_parameters, model.grass_parameters, model.rip_parameters, model.slow_parameters, model.Total_Elevationbands, model.Elevation_Percentage)
    end 
    
    #Update model state
    model.Discharge = Discharge
    model.Snow_Extend = Snow_Extend
    model.Waterbalance = Waterbalance
    model.bare_storage = bare_storage::Storages
    model.forest_storage = forest_storage::Storages
    model.grass_storage = grass_storage::Storages
    model.rip_storage = rip_storage::Storages
    model.Slowstorage = Slowstorage::Float64
    
    #update current date
    model.Current_Date = model.Date[t+1]

end 

# Updates model until indicated time step
function update_until(model::HBVmountain_model, time::String)
    update_time = length(model.Date[model.Date .<= DateTime(time)])
    tmax::Int128 = length(model.Precipitation[1:update_time])
    for t in 1:tmax
        update(model)
    end
end

# Returns empty model
function finalize(model) 
    finalize_HBVmountainmodel(model)
end

function get_component_name(model)
    return string(typeof(model))
end

function get_input_item_count(model)
    return length(fieldnames(typeof(model)))
end

function get_output_item_count(model)
    return length(fieldnames(typeof(model)))
end

function get_input_var_names(model)
    return string(fieldnames(typeof(model)))
end

function get_output_var_names(model)
    return string(fieldnames(typeof(model)))
end

function get_var_type(model, name::String)
    return String(typeof(getfield(model, Symbol(name))))
end

function get_var_units(model, name::String)
    return getfield(model.Units, Symbol(name))
end

function get_var_itemsize(model, name::String)
    return Int(sizeof(getfield(model, Symbol(name)))/ length(getfield(model, Symbol(name))))
end

function get_var_nbytes(model, name::String)
    return sizeof(getfield(model, Symbol(name)))
end

function get_current_time(model)
    return model.Current_Date
end

function get_start_time(model)
    return model.Date[1]
end

function get_end_time(model)
    return last(model.Date)
end

function get_time_units(model)
    return dump(model.Date[2] - model.Date[1], maxdepth=0)
end

function get_time_step(model)
    return (model.Date[2] - model.Date[1])::Float64
end

function get_value(model, name, dest=nothing)
    if dest != nothing 
        empty!(dest)
        push!(dest, getfield(model, Symbol(name)))
    else   
        return getfield(model, Symbol(name))
    end 
end

# Not sure if this is correct
function get_value_ptr(model, name)
    return getfield(model, Symbol(name))
end

function get_value_at_indices(model, name, dest, inds)
    getfield(model, Symbol(name))[inds]
end

function set_value(model, name::String, value)
    setfield!(model, Symbol(name), value)
end

function set_value_at_indices(model, name, inds, value)
    model.name[inds] = value
end

function get_var_grid(model, name)
    error("BMI Function not implemented"::AbstractString) 
end

function get_var_location(model, name)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_rank(model, grid)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_size(model, grid)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_type(model, grid)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_shape(model, grid)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_spacing(model, grid, spacing)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_origin(model, grid, origin)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_x(model, grid, x)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_y(model, grid, y)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_z(model, grid, z)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_node_count(model, grid)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_edge_count(model, grid)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_face_count(model, grid)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_edge_nodes(model, grid, edge_nodes)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_face_edges(model, grid, face_edges)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_face_nodes(model, grid, face_nodes)
    error("BMI Function not implemented"::AbstractString) 
end

function get_grid_nodes_per_face(model, grid, nodes_per_face)
    error("BMI Function not implemented"::AbstractString) 
end















