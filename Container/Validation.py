from Calibration import *


####### Forcing and observation data ##########################
def run_validation(calibration_results, path_to_observation, forcing, path_to_shapefile, path_to_dem, path_to_nlcd, end_year_calibration, end_year_validation):
    """
    Function for validation of the HBV-mountain model. This function also returns the streamflow simulations of the total calibration and validation period, together with the simulated evaporation.
    
    Parameters
    ------------
    calibration_results: Dataframe generated using the calibration function. Contains the parameter set and corresponding objective scores.
    path_to_observation: str. Path to GRDC streamflow observation data. This must be a CSV file containing column called 'streamflow' containing the streamflow data.
    forcing: netCDF4.Dataset containing precipitation and temperature data
    path_to_shapefile: str. Path to catchment shapefile (must be WGS84)
    path_to_dem: str. Path to DEM raster file (must be WGS84)
    path_to_nlcd: str. Path to NLCD raster file (must be WGS84)
    end_year_calibration: Int. Last year of the calibration period
    end_year_validation: Int. Last year of validation period
    
    
    Returns Dataframe with parameter set(s) and corresponding objective scores in the validation period as columns, and the results of the different validation runs as rows. Also returns the streamflow simulations of the total calibration and validation period, together with the simulated evaporation.
    """
    
    
    ob_list = []
    params_list = []
    sim_list = []
    df_evap = pd.DataFrame(index=generate_forcing_from_NETCDF(forcing).prec.index)
    glacier = []
    df = pd.DataFrame(index=generate_forcing_from_NETCDF(forcing).prec.index)
    for i in range(len(calibration_results)):
        parameters = calibration_results.iloc[i, -20:]
        model = BMI_HBVmountain(forcing, path_to_shapefile, path_to_dem, path_to_nlcd)
        config_file = model.setup(forcing_netcdf=forcing, bare_parameters=  Parameters(parameters.beta_Bare, parameters.Ce, 0, 0, parameters.Kf,
                                                                         parameters.Meltfactor, parameters.Mm, parameters.Ratio_Pref, parameters.Soilstoragecapacity_Bare, parameters.Temp_Thresh),
                                            forest_parameters=Parameters(parameters.beta_Forest, parameters.Ce, 0, parameters.Interceptioncapacity_Forest, parameters.Kf,
                                                                         parameters.Meltfactor, parameters.Mm, parameters.Ratio_Pref, parameters.Soilstoragecapacity_Forest, parameters.Temp_Thresh),
                                            grass_parameters= Parameters(parameters.beta_Grass, parameters.Ce, 0, parameters.Interceptioncapacity_Grass, parameters.Kf,
                                                                         parameters.Meltfactor, parameters.Mm, parameters.Ratio_Pref, parameters.Soilstoragecapacity_Grass, parameters.Temp_Thresh),
                                            rip_parameters=   Parameters(parameters.beta_Rip, parameters.Ce, 0, parameters.Interceptioncapacity_Rip, parameters.Kf_Rip,
                                                                         parameters.Meltfactor, parameters.Mm, parameters.Ratio_Pref, parameters.Soilstoragecapacity_Rip, parameters.Temp_Thresh),
                                            slow_parameters=  Slow_Paramters(parameters.Ks, parameters.Ratio_Riparian))
        model.initialize(config_file)



    ########### Catchment dependent settings ######################


        Evaporation = []
        Discharge = []
        timestamp = []
        while (model.get_value_ptr('Current_Date') < (datetime.date(end_year_validation, 12, 31))):  
            model.update()
            timestamp.append(model.get_value_ptr('Current_Date'))
            Discharge.append(model.get_value_ptr('Discharge'))
            Evaporation.append(model.get_value_ptr('Total_Evaporation'))
            glacier.append(model.get_value_ptr('Glacier'))

        simulated_discharge = simulated_discharge_df =  pd.DataFrame(
                {'streamflow': Discharge},
                index=pd.to_datetime(timestamp)
            )
        simulated_discharge_df = pd.DataFrame(
            {f'simulation_{i+1}': Discharge},
            index=pd.to_datetime(timestamp)
        )
        simulated_evaporation_df = pd.DataFrame(
            {f'simulation_{i+1}': Evaporation},
            index=pd.to_datetime(timestamp)
        )
        df = pd.merge(df, simulated_discharge_df, left_index=True, right_index=True)
        df_evap = pd.merge(df_evap, simulated_evaporation_df, left_index=True, right_index=True)
        
        model.finalize()
        
        # Validation 
        precipitation = generate_forcing_from_NETCDF(forcing).prec
        area = gpd.read_file(path_to_shapefile).area_hys.values[0] #km2 
        observation = pd.read_csv(path_to_observation, index_col=0).streamflow / (area * 1e6) * 1000 *86400 
        mask = (observation.index >= str(simulated_discharge.index[0])) & (observation.index <= str(simulated_discharge.index[-1]))
        observation = observation.loc[mask]
        observation.index = pd.to_datetime(observation.index)
        precipitation.index = pd.to_datetime(precipitation.index)
        objective_function = multi_objective(simulated_discharge.loc[simulated_discharge.index.year > end_year_calibration].streamflow, 
                                             observation.loc[observation.index.year > end_year_calibration], precipitation)
        ob_list.append(objective_function)
        params_list.append(parameters)    

        paramset = pd.DataFrame(columns=['ED_val', 'NSE_val', 'logNSE_val', 'NSEfdc_val', 'NSErunoff_val', 'Temp_Thresh', 'Meltfactor', 'Mm', 'Ratio_Pref', 'Kf', 'Ks', 'Ce', 'Soilstoragecapacity_Bare', 'beta_Bare', 'Interceptioncapacity_Forest', 'Soilstoragecapacity_Forest', 'beta_Forest', 'Interceptioncapacity_Grass', 'Soilstoragecapacity_Grass', 'beta_Grass', 'Interceptioncapacity_Rip', 'Soilstoragecapacity_Rip', 'beta_Rip', 'Kf_Rip', 'Ratio_Riparian'])
        
        
    for i in range(len(ob_list)):
        paramset.loc[i] = [ob_list[i][0], ob_list[i][1], ob_list[i][2], ob_list[i][3], ob_list[i][4], params_list[i][0], params_list[i][1], params_list[i][2], params_list[i][3], params_list[i][4], params_list[i][5], params_list[i][6], params_list[i][7], params_list[i][8], params_list[i][9], params_list[i][10], params_list[i][11], params_list[i][12], params_list[i][13], params_list[i][14], params_list[i][15], params_list[i][16], params_list[i][17], params_list[i][18], params_list[i][19]]    
    paramset.ED_val = 1 - paramset.ED_val

    return paramset, df, df_evap
