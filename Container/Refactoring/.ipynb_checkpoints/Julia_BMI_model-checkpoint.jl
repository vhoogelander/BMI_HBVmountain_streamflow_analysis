using DataFrames
using CSV
using DocStringExtensions
using Dates
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
function setup(Discharge, Total_Evaporation, Snow_Extend, bare_storage,forest_storage, grass_storage, rip_storage,
        Slowstorage, Waterbalance, Glacier, Area, bare_parameters, forest_parameters, grass_parameters, rip_parameters, slow_parameters,
        Elevation, Total_Elevationbands, Precipitation_gradient, Elevation_Percentage, bare_input, forest_input, grass_input, rip_input, Precipitation, Temperature, ETP, Date, Current_Date, Sunhours, Units)
        """
        Returns a configuration file for the setup of BMI HBV-mountain model.
        """
    config_file = HBVmountain_model(Discharge, Total_Evaporation, Snow_Extend, bare_storage,forest_storage, grass_storage, rip_storage,
        Slowstorage, Waterbalance, Glacier, Area, bare_parameters, forest_parameters, grass_parameters, rip_parameters, slow_parameters,
        Elevation, Total_Elevationbands, Precipitation_gradient, Elevation_Percentage, bare_input, forest_input, grass_input, rip_input, Precipitation, Temperature, ETP, Date, Current_Date, Sunhours, Units);
    return config_file
end


function initialize(model, config_file = nothing)
    """
    Initializes BMI HBV-mountain model using a configuration file
    """
    if config_file != nothing
        replace_HBVmountainmodel(model, config_file) 
    else
        error("must give config_file"::AbstractString)
    end
end

# Updates model for one time step, uses the function in Refactoring/run_model 
function update(model::HBVmountain_model)
    """
    Updates the BMI HBV-mountain model for one timestep
    """    
    
    
    
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
        Discharge, Total_Evaporation, Snow_Extend, Waterbalance, bare_storage, forest_storage, grass_storage, rip_storage, Slowstorage =  run_model_glacier(model.Area, model.ETP[t], model.Glacier, Precipitation, Temperature,
                        model.bare_input, model.forest_input, model.grass_input, model.rip_input,
                        model.bare_storage, model.forest_storage, model.grass_storage, model.rip_storage, model.Slowstorage,
                        model.bare_parameters, model.forest_parameters, model.grass_parameters, model.rip_parameters, model.slow_parameters, model.Total_Elevationbands, model.Elevation_Percentage)
    end 
    
    #Update model state
    model.Discharge = Discharge
    model.Total_Evaporation = Total_Evaporation
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
    """
    Updates the BMI HBV-mountain model for until given time. Time must be a string.
    """
    update_time = length(model.Date[model.Date .<= DateTime(time)])
    tmax::Int128 = length(model.Precipitation[1:update_time])
    for t in 1:tmax
        update(model)
    end
end

# Returns empty model
function finalize(model) 
    """
    Perform tear-down tasks for the model.
    """
    finalize_HBVmountainmodel(model)
end

function get_component_name(model)
    """
    Get name of the component
    """
    return string(typeof(model))
end

function get_input_item_count(model)
    """
    Count of a model's input variables.
    
    Returns
    -------
    int
    """
    return length(fieldnames(typeof(model)))
end

function get_output_item_count(model)
    """
    Count of a model's output variables.
    
    Returns
    -------
    int
    """
    return length(fieldnames(typeof(model)))
end

function get_input_var_names(model)
    """
    List of a model's input variables.
    
    Returns
    -------
    list of str
    """

    return string(fieldnames(typeof(model)))
end

function get_output_var_names(model)
    """
    List of a model's output variables.
    
    Returns
    -------
    list of str
    """
    return string(fieldnames(typeof(model)))
end

function get_var_type(model, name::String)
    """
    Get data type of the given variable.
    Parameters
    ----------
    name : str
        An input or output variable name
    Returns
    -------
    str
    """
    return String(typeof(getfield(model, Symbol(name))))
end

function get_var_units(model, name::String)
    """
    Get units of the given variable
    Parameters
    ----------
    name : str
        An input or output variable name
    Returns
    -------
    str
    """
    return getfield(model.Units, Symbol(name))
end

function get_var_itemsize(model, name::String)
    """
    Get memory use for each array element in bytes.
    Parameters
    ----------
    name : str
        An input or output variable name
    Returns
    -------
    int
    """
    return Int(sizeof(getfield(model, Symbol(name)))/ length(getfield(model, Symbol(name))))
end

function get_var_nbytes(model, name::String)
    """
    Get size, in bytes, of the given variable.
    Parameters
    ----------
    name : str
        An input or output variable name
    Returns
    -------
    int
    """
    return sizeof(getfield(model, Symbol(name)))
end

function get_current_time(model)
    """
    Current time of the model.

    Returns
    -------
    datetime
    """
    return model.Current_Date
end

function get_start_time(model)
    """
    Start time of the model.

    Returns
    -------
    datetime
    """
    return model.Date[1]
end

function get_end_time(model)
    """
    End time of the model.

    Returns
    -------
    datetime
    """
    return last(model.Date)
end

function get_time_units(model)
    """
    Time units of the model.

    Returns
    -------
    String
    """
    return string(typeof((model.Date[2] - model.Date[1])))
end

function get_time_step(model)
    """
    Current time step of the model.

    Returns
    -------
    float
    """
    return (model.Date[2] - model.Date[1]) / (model.Date[2] - model.Date[1])
end

function get_value(model, name, dest=nothing)
    """
    Get a copy of values of the given variable.

    This is a getter for the model, used to access the model's
    current state. It returns a *copy* of a model variable, with
    the return type, size and rank dependent on the variable.

    Parameters
    ----------
    name : str
        An input or output variable name
    dest : ndarray
        An array into which to place the values.

    Returns
    -------
    array
    """
    if dest != nothing 
        empty!(dest)
        push!(dest, getfield(model, Symbol(name)))
    else   
        return getfield(model, Symbol(name))
    end 
end


function get_value_ptr(model, name)
    """
    Get a reference to values of the given variable.

    This is a getter for the model, used to access the model's
    current state. It returns a reference to a model variable,
    with the return type, size and rank dependent on the variable.

    Parameters
    ----------
    name : str
        An input or output variable name

    Returns
    -------
    array_like
    """
    return getfield(model, Symbol(name))
end

function get_value_at_indices(model, name, dest, inds)
    """
    Get values at particular indices.

    Parameters
    ----------
    name : str
        An input or output variable name
    dest : ndarray
        An array into which to place the values.
    indices : array_like
        The indices into the variable array.

    Returns
    -------
    array_like
        Value of the model variable at the given location.
    """
    getfield(model, Symbol(name))[inds]
end

function set_value(model, name::String, value)
    """
    Specify a new value for a model variable.

    This is the setter for the model, used to change the model's
    current state. It accepts, through *src*, a new value for a
    model variable, with the type, size and rank of *src*
    dependent on the variable.

    Parameters
    ----------
    name : str
        An input or output variable name
    value : array_like
        The new value for the specified variable.
    """
    setfield!(model, Symbol(name), value)
end

function set_value_at_indices(model, name, inds, value)
    """
    Specify a new value for a model variable at particular indices.

    This is the setter for the model, used to change the model's
    current state. It accepts, through *src*, a new value for a
    model variable, with the type, size and rank of *src*
    dependent on the variable.

    Parameters
    ----------
    name : str
        An input or output variable name
    indices : array_like
        The indices into the variable array.
    value : array_like
        The new value for the specified variable.
    """
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















