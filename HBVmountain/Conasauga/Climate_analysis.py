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

        model.set_value('Elevation', Main.Elevations(274.25, 187.0, 1284.0, 421.17, 421.17))
        model.set_value('Glacier', [0.0, 0.0, 0.0, 0.0])
        model.set_value('Sunhours', [10.15, 11.0, 11.93, 13.00, 13.9, 14.38, 14.18, 13.42, 12.38, 11.35, 10.42, 9.92])
        model.set_value('bare_input', Main.HRU_Input([1.0,0.0,0.0,0.0], 0.055, np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('forest_input', Main.HRU_Input([0.58,0.22,0.16,0.04], 0.7,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('grass_input', Main.HRU_Input([0.97, 0.02, 0.01, 0.0], 0.214,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('rip_input', Main.HRU_Input([0.99,0.01,0.0,0.0], 0.001,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('Total_Elevationbands', 4)
        model.set_value('Elevation_Percentage', [0.68,0.168,0.125,0.027])


        Discharge = []
        timestamp = []
        while (model.get_value_ptr('Current_Date') < (model.get_value_ptr('Date'))[-1]):  
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
        model.set_value('Elevation', Main.Elevations(274.25, 187.0, 1284.0, 421.17, 421.17))
        model.set_value('Glacier', [0.0, 0.0, 0.0, 0.0])
        model.set_value('Sunhours', [10.15, 11.0, 11.93, 13.00, 13.9, 14.38, 14.18, 13.42, 12.38, 11.35, 10.42, 9.92])
        model.set_value('bare_input', Main.HRU_Input([1.0,0.0,0.0,0.0], 0.055, np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('forest_input', Main.HRU_Input([0.58,0.22,0.16,0.04], 0.7,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('grass_input', Main.HRU_Input([0.97, 0.02, 0.01, 0.0], 0.214,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('rip_input', Main.HRU_Input([0.99,0.01,0.0,0.0], 0.001,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
        model.set_value('Total_Elevationbands', 4)
        model.set_value('Elevation_Percentage', [0.68,0.168,0.125,0.027])



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