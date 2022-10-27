from Validation import *


def linear_scaling_correction(str_path_calibration_forcing, str_path_input_hist_forcing, str_path_input_forcing):
    calibration_forcing = nc.Dataset(str_path_calibration_forcing)
    input_hist_forcing = nc.Dataset(str_path_input_hist_forcing)
    input_forcing = nc.Dataset(str_path_input_forcing)

    calibration_forcing  = generate_forcing_from_NETCDF(calibration_forcing)
    calibration_forcing.index = pd.to_datetime(calibration_forcing.index)
    calibration_forcing = calibration_forcing.loc[calibration_forcing.index.year <= 2005]

    input_hist_forcing = generate_forcing_from_NETCDF(input_hist_forcing)
    input_hist_forcing.index = pd.to_datetime(input_hist_forcing.index)
    input_hist_forcing = input_hist_forcing.loc[(input_hist_forcing.index.year >= calibration_forcing.index.year[0]) & (input_hist_forcing.index.year <= calibration_forcing.index.year[-1])]

    daily_mean_prec_input = input_hist_forcing.groupby(input_hist_forcing.index.strftime("%m")).mean() #.prec.transform('mean')
    daily_mean_prec_calibration = calibration_forcing.groupby(calibration_forcing.index.strftime("%m")).mean() #.prec.transform('mean')
    
    prec_correction_factor = daily_mean_prec_calibration.prec / daily_mean_prec_input.prec
    prec_correction_factor.index = np.arange(1,13)
    
    temp_correction_factor = daily_mean_prec_calibration.temp - daily_mean_prec_input.temp
    temp_correction_factor.index = np.arange(1,13)

    
    input_forcing  = generate_forcing_from_NETCDF(input_forcing)
    input_forcing = nc.Dataset(str_path_input_forcing)
    input_forcing  = generate_forcing_from_NETCDF(input_forcing)
    input_forcing.index = pd.to_datetime(input_forcing.index)

    for m in range(len(prec_correction_factor)):
        input_forcing.loc[input_forcing.index.month == m+1, 'prec'] = input_forcing.loc[input_forcing.index.month == m+1].prec * prec_correction_factor.values[m]
        input_forcing.loc[input_forcing.index.month == m+1, 'temp'] = input_forcing.loc[input_forcing.index.month == m+1].temp + temp_correction_factor.values[m]


    
    return input_forcing
    
def run_climate_simulations(calval_results, str_path_forcing):
    forcing = nc.Dataset(str_path_forcing) #Catchment dependent
    ob_list = []
    params_list = []
    sim_list = []
    
    df = pd.DataFrame(index=generate_forcing_from_NETCDF(forcing).prec.index)
    for i in range(len(calval_results)):
        parameters = calval_results.iloc[i, -20:]
        model = BMI_HBVmountain(forcing_netcdf=forcing)
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

        model.set_value('Elevation', Main.Elevations(386.75, 1245.0, 2792.0, 1902.188, 1902.188))
        model.set_value('Glacier', [0.0, 0.0, 0.0, 0.0])
        model.set_value('Sunhours', [10.2 ,11.0 ,12.0 ,13.1 ,14.0 ,14.4, 14.1 ,13.3, 12.3 ,11.3, 10.4, 9.9])
        model.set_value('bare_input', Main.HRU_Input([0.277, 0.161, 0.127, 0.435], 0.023, np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('forest_input', Main.HRU_Input([0.113, 0.466, 0.304, 0.117], 0.623,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('grass_input', Main.HRU_Input([0.491, 0.261, 0.192, 0.056], 0.354,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('rip_input', Main.HRU_Input([0.584, 0.37, 0.033, 0.013], 0.01,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('Total_Elevationbands', 4)
        model.set_value('Elevation_Percentage', [0.254, 0.385, 0.269, 0.092])



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

def run_climate_simulations_biascorrected(calval_results, str_path_input_forcing, str_path_calibration_forcing, str_path_input_hist_forcing):
    forcing = nc.Dataset(str_path_input_forcing)
    forcing_corrected = linear_scaling_correction(str_path_calibration_forcing, str_path_input_hist_forcing, str_path_input_forcing) 
    ob_list = []
    params_list = []
    sim_list = []
    
    df = pd.DataFrame(index=generate_forcing_from_NETCDF(forcing).prec.index)
    for i in range(len(calval_results)):
        parameters = calval_results.iloc[i, -20:]
        model = BMI_HBVmountain(forcing_netcdf=forcing)
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
        model.set_value('Elevation', Main.Elevations(386.75, 1245.0, 2792.0, 1902.188, 1902.188))
        model.set_value('Glacier', [0.0, 0.0, 0.0, 0.0])
        model.set_value('Sunhours', [10.2 ,11.0 ,12.0 ,13.1 ,14.0 ,14.4, 14.1 ,13.3, 12.3 ,11.3, 10.4, 9.9]) #Los Angeles
        model.set_value('bare_input', Main.HRU_Input([0.277, 0.161, 0.127, 0.435], 0.023, np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('forest_input', Main.HRU_Input([0.113, 0.466, 0.304, 0.117], 0.623,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('grass_input', Main.HRU_Input([0.491, 0.261, 0.192, 0.056], 0.354,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('rip_input', Main.HRU_Input([0.584, 0.37, 0.033, 0.013], 0.01,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('Total_Elevationbands', 4)
        model.set_value('Elevation_Percentage', [0.254, 0.385, 0.269, 0.092])



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