from Validation import *

####### Forcing and observation data ##########################
def run_climate_simulations(calval_results, str_path_input_forcing, path_to_shapefile, path_to_dem, path_to_nlcd):
    """
    Function for running streamflow simulations
    
    Parameters
    --------------
    calval_results: Dataframe containing parameter sets from calibration/validation process
    str_path_input_forcing: str. Historical data of climate model
    path_to_shapefile: str. Path to catchment shapefile (must be WGS84)
    path_to_dem: str. Path to DEM raster file (must be WGS84)
    path_to_nlcd: str. Path to NLCD raster file (must be WGS84)
    
    
    Returns Dataframe containing streamflow simulations, with n columns of the simulations of n different parameter sets
    """
    forcing = nc.Dataset(str_path_input_forcing)
    
    df = pd.DataFrame(index=generate_forcing_from_NETCDF(forcing).prec.index)
    for i in range(len(calval_results)):
        parameters = calval_results.iloc[i]
        model = BMI_HBVmountain(forcing, path_to_shapefile, path_to_dem, path_to_nlcd)
        config_file = model.setup(bare_parameters=  Parameters(parameters.beta_Bare, parameters.Ce, 0, 0, parameters.Kf,
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

        Discharge = []
        timestamp = []
        while (model.get_value_ptr('Current_Date') < (model.get_value_ptr('Date')[-1])):  
            model.update()
            timestamp.append(model.get_value_ptr('Current_Date'))
            Discharge.append(model.get_value_ptr('Discharge'))

        simulated_discharge_df = pd.DataFrame(
            {f'simulation_{i+1}': Discharge},
            index=pd.to_datetime(timestamp)
        )
        df = pd.merge(df, simulated_discharge_df, left_index=True, right_index=True)
        
        model.finalize()
        
    
    return df


def run_climate_simulations_biascorrected(calval_results, str_path_input_forcing, str_path_calibration_forcing, str_path_input_hist_forcing, path_to_shapefile, path_to_dem, path_to_nlcd):
    """
    Function for running streamflow simulations, where the forcing is bias corrected with the calibration forcing. Function tested and used for CMIP6 forcing data and ERA5 calibration data.
    
    Parameters
    --------------
    calval_results: Dataframe containing parameter sets from calibration/validation process
    str_path_input_forcing: str. Historical data of climate model
    str_path_calibration_forcing: str. Data used for model calibration
    str_path_input_hist_forcing:  str. Climate data used for streamflow simulations
    path_to_shapefile: str. Path to catchment shapefile (must be WGS84)
    path_to_dem: str. Path to DEM raster file (must be WGS84)
    path_to_nlcd: str. Path to NLCD raster file (must be WGS84)
    
    
    Returns Dataframe containing streamflow simulations, with n columns of the simulations of n different parameter sets
    """
    
    forcing = nc.Dataset(str_path_input_forcing)
    forcing_corrected = linear_scaling_correction(str_path_calibration_forcing, str_path_input_hist_forcing, str_path_input_forcing) 
    
    df = pd.DataFrame(index=generate_forcing_from_NETCDF(forcing).prec.index)
    for i in range(len(calval_results)):
        parameters = calval_results.iloc[i]
        model = BMI_HBVmountain(forcing, path_to_shapefile, path_to_dem, path_to_nlcd)
        config_file = model.setup(bare_parameters=  Parameters(parameters.beta_Bare, parameters.Ce, 0, 0, parameters.Kf,
                                                                         parameters.Meltfactor, parameters.Mm, parameters.Ratio_Pref, parameters.Soilstoragecapacity_Bare, parameters.Temp_Thresh),
                                            forest_parameters=Parameters(parameters.beta_Forest, parameters.Ce, 0, parameters.Interceptioncapacity_Forest, parameters.Kf,
                                                                         parameters.Meltfactor, parameters.Mm, parameters.Ratio_Pref, parameters.Soilstoragecapacity_Forest, parameters.Temp_Thresh),
                                            grass_parameters= Parameters(parameters.beta_Grass, parameters.Ce, 0, parameters.Interceptioncapacity_Grass, parameters.Kf,
                                                                         parameters.Meltfactor, parameters.Mm, parameters.Ratio_Pref, parameters.Soilstoragecapacity_Grass, parameters.Temp_Thresh),
                                            rip_parameters=   Parameters(parameters.beta_Rip, parameters.Ce, 0, parameters.Interceptioncapacity_Rip, parameters.Kf_Rip,
                                                                         parameters.Meltfactor, parameters.Mm, parameters.Ratio_Pref, parameters.Soilstoragecapacity_Rip, parameters.Temp_Thresh),
                                            slow_parameters=  Slow_Paramters(parameters.Ks, parameters.Ratio_Riparian))
        model.initialize(config_file)


        model.set_value('Temperature', (np.float64(forcing_corrected['temp'].values)).reshape(len(forcing_corrected),1))
        model.set_value('Precipitation', (np.float64(forcing_corrected['prec'].values)).reshape(len(forcing_corrected),1))
    ########### Catchment dependent settings ######################

        Discharge = []
        timestamp = []
        while (model.get_value_ptr('Current_Date') < (model.get_value_ptr('Date')[-1])):  
            model.update()
            timestamp.append(model.get_value_ptr('Current_Date'))
            Discharge.append(model.get_value_ptr('Discharge'))

        simulated_discharge_df = pd.DataFrame(
            {f'simulation_{i+1}': Discharge},
            index=pd.to_datetime(timestamp)
        )
        df = pd.merge(df, simulated_discharge_df, left_index=True, right_index=True)
        
        model.finalize()
    
    return df

def create_monthly_boxplots(simulations_hist, simulations_ssp245, simulations_ssp585):
    fig, axarr = plt.subplots(figsize=(12,4))

    hist = simulations_hist.groupby(simulations_hist.index.strftime("%y-%m")).sum()
    hist['mean'] = hist.mean(axis=1)
    hist['month'] = hist.index.str[3:]
    hist.boxplot(by='month', column='mean', ax=axarr, positions=np.array(range(12))*3.0-0.8, sym='', widths=0.6, color='k')

    ssp245 = simulations_ssp245.groupby(simulations_ssp245.index.strftime("%y-%m")).sum()
    ssp245['mean'] = ssp245.mean(axis=1)
    ssp245['month'] = ssp245.index.str[3:]
    ssp245.boxplot(by='month', column='mean', ax=axarr, sym='', positions=np.array(range(12))*3.0, widths=0.6, color='b')

    ssp585 = simulations_ssp585.groupby(simulations_ssp585.index.strftime("%y-%m")).sum()
    ssp585['mean'] = ssp585.mean(axis=1)
    ssp585['month'] = ssp585.index.str[3:]
    ssp585.boxplot(by='month', column='mean', ax=axarr, positions=np.array(range(12))*3.0+0.8, sym='', widths=0.6, color='r')
    
    ticks = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    fig.suptitle('')
    plt.xticks(range(0, 12 * 3, 3), ticks)