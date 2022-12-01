from BMI_HBVmountain_Python import *
import os
import sys
import netCDF4 as nc
####### Forcing and observation data ##########################

def NSE(Qmodelled, Qobserved):
    QobservedAverage = np.ones(len(Qobserved)) * np.mean(Qobserved)
    Nominator = np.sum((Qmodelled - Qobserved)**2)
    Denominator = np.sum((Qobserved - QobservedAverage)**2)
    NashSutcliffe = 1 - (Nominator / Denominator)
    return NashSutcliffe

def logNSE(Qmodelled, Qobserved):
    QobservedAverage = np.ones(len(Qobserved)) * np.mean(Qobserved) # average as array
    Nominator = np.sum((np.log(Qobserved)-np.log(Qmodelled))**2)
    Denominator = np.sum((np.log(Qobserved) - np.log(QobservedAverage))**2)
    NashSutcliffelog = 1 - (Nominator / Denominator)
    return NashSutcliffelog

def flowdurationcurve(Q):
    SortedQ = np.sort(Q)[::-1]
    rank = np.linspace(1, len(SortedQ), len(SortedQ))
    ExcProb = rank / (len(SortedQ) + 1)
    return SortedQ, ExcProb

# def autocorrelation(Q, Timelag)
#     Qshifted = Q[1 + Timelag: end]
#     Q = Q[1 : end - Timelag]
#     Q_average = ones(length(Q)) * mean(Q)
#     Nominator = sum((Q -  Q_average) .* (Qshifted - Q_average))
#     Denominator = sum((Q - Q_average).^2)
#     AC = Nominator / Denominator
#     return AC::Float64

def yearlyrunoff(Precipitation, Discharge):
    annual_prec = Precipitation.groupby(pd.PeriodIndex(Precipitation.index, freq="y")).sum()
    annual_Discharge = Discharge.groupby(pd.PeriodIndex(Discharge.index, freq="y")).sum()
    mask = (annual_prec.index >= annual_Discharge.index[0]) & (annual_prec.index <= annual_Discharge.index[-1])
    annual_prec = annual_prec.loc[mask]
    annual_runoff_coefficient = annual_Discharge / annual_prec
    return annual_runoff_coefficient

def multi_objective(Qmodelled, Qobserved, Precipitation):
    nse = NSE(Qmodelled.values, Qobserved.values)
    lognse = logNSE(Qmodelled.values, Qobserved.values)

    FDCobserved = flowdurationcurve(np.log(Qobserved))
    FDCmodelled = flowdurationcurve(np.log(Qmodelled))
    nseFDC = NSE(FDCmodelled[0], FDCobserved[0])

    runoff_observed = yearlyrunoff(Precipitation, Qobserved)
    runoff_modelled = yearlyrunoff(Precipitation, Qmodelled)
    nse_runoff = NSE(runoff_modelled.values, runoff_observed.values)

    objfunctions = [nse, lognse, nseFDC, nse_runoff]
    Sum = 0
    for obj in objfunctions:
        Sum += (1 - obj)**2
    ED = (Sum / len(objfunctions))**0.5

    return ED, nse, lognse, nseFDC, nse_runoff

def parameter_conversion(directory, nsmallest=None):
    ## Returns a list with Julia parameter objects 
    bare_parametersets = pd.read_csv(os.path.join(directory, 'bare_paramsets.csv')).sort_values(by=['ED'])
    if nsmallest != None:
        bare_parametersets = bare_parametersets.loc[bare_parametersets['ED'].isin(bare_parametersets.ED.nsmallest(nsmallest))]
    forest_parametersets = pd.read_csv(os.path.join(directory, 'forest_paramsets.csv')).loc[bare_parametersets.index].sort_values(by=['ED'])
    grass_parametersets = pd.read_csv(os.path.join(directory, 'grass_paramsets.csv')).loc[bare_parametersets.index].sort_values(by=['ED'])
    rip_parametersets = pd.read_csv(os.path.join(directory, 'rip_paramsets.csv')).loc[bare_parametersets.index].sort_values(by=['ED'])
    slow_parametersets = pd.read_csv(os.path.join(directory, 'slow_paramsets.csv')).loc[bare_parametersets.index].sort_values(by=['ED'])
    
    params = [bare_parametersets, forest_parametersets, grass_parametersets, rip_parametersets]
    param_list = []
    for i in range(len(bare_parametersets)):
        bare = Parameters(bare_parametersets.beta.iloc[i], bare_parametersets.Ce.iloc[i], bare_parametersets.Drainagecapacity.iloc[i], 
                                             bare_parametersets.Interceptionstoragecapacity.iloc[i], bare_parametersets.Kf.iloc[i], 
                                             bare_parametersets.Meltfactor.iloc[i], bare_parametersets.Mm.iloc[i], bare_parametersets.Ratio_Pref.iloc[i], 
                                            bare_parametersets.Soilstoragecapacity.iloc[i], bare_parametersets.Temp_Thresh.iloc[i])
        forest = Parameters(forest_parametersets.beta.iloc[i], forest_parametersets.Ce.iloc[i], forest_parametersets.Drainagecapacity.iloc[i], 
                                            forest_parametersets.Interceptionstoragecapacity.iloc[i], forest_parametersets.Kf.iloc[i], 
                                            forest_parametersets.Meltfactor.iloc[i], forest_parametersets.Mm.iloc[i], forest_parametersets.Ratio_Pref.iloc[i], 
                                            forest_parametersets.Soilstoragecapacity.iloc[i], forest_parametersets.Temp_Thresh.iloc[i])
        grass = Parameters(grass_parametersets.beta.iloc[i], grass_parametersets.Ce.iloc[i], grass_parametersets.Drainagecapacity.iloc[i], 
                                             grass_parametersets.Interceptionstoragecapacity.iloc[i], grass_parametersets.Kf.iloc[i], 
                                             grass_parametersets.Meltfactor.iloc[i], grass_parametersets.Mm.iloc[i], grass_parametersets.Ratio_Pref.iloc[i], 
                                            grass_parametersets.Soilstoragecapacity.iloc[i], grass_parametersets.Temp_Thresh.iloc[i])
        rip = Parameters(rip_parametersets.beta.iloc[i], rip_parametersets.Ce.iloc[i], rip_parametersets.Drainagecapacity.iloc[i], 
                                             rip_parametersets.Interceptionstoragecapacity.iloc[i], rip_parametersets.Kf.iloc[i], 
                                             rip_parametersets.Meltfactor.iloc[i], rip_parametersets.Mm.iloc[i], rip_parametersets.Ratio_Pref.iloc[i], 
                                            rip_parametersets.Soilstoragecapacity.iloc[i], rip_parametersets.Temp_Thresh.iloc[i])
        slow = Slow_Paramters(slow_parametersets.Ks.iloc[i], slow_parametersets.Ratio_Riparian.iloc[i])
    
        param_list.append([bare, forest, grass, rip, slow])
    return param_list

def climate_simulations(prec, temp, param_list):
    # ## Setting up the model
    forcing = pd.DataFrame({'prec':prec, 'temp':temp}, index=prec.index)

    df = pd.DataFrame(index=forcing.index)
    for i in range(len(param_list)):

        model = BMI_HBVmountain()
        config_file = model.setup()
        model.initialize(config_file)


        # ### Parameters
        model.set_value('bare_parameters', param_list[i][0])
        model.set_value('forest_parameters', param_list[i][1])
        model.set_value('grass_parameters', param_list[i][2])
        model.set_value('rip_parameters', param_list[i][3])
        model.set_value('slow_parameters', param_list[i][4])


        #Historical forcing
        model.set_value('Temperature', (forcing['temp'].values).reshape(len(forcing),1))
        model.set_value('Precipitation', (forcing['prec'].values).reshape(len(forcing),1))

        model.set_value('Date', list(forcing.index.values))
        model.set_value('Current_Date', forcing.index.values[0])



        # ### Initial settings
        model.set_value('Elevation', Main.Elevations(270, 1250, 2600, 1250, 1250))
        model.set_value('Glacier', [0.0, 0.0, 0.0, 0.0])
        model.set_value('Sunhours', [10.18, 10.90, 12.0, 13.10, 14.0, 14.45, 14.08, 13.31, 12.24, 11.25, 10.31, 9.85])
        model.set_value('bare_input', Main.HRU_Input([0,0,0,0], 0, np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.2, np.zeros(4), 0, 0.0))
        model.set_value('forest_input', Main.HRU_Input([0.0,0.20,0.73,0.07], 0.64,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.2, np.zeros(4), 0, 0.0))
        model.set_value('grass_input', Main.HRU_Input([0.1,0.8,0.1,0.0], 0.35,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.2, np.zeros(4), 0, 0.0))
        model.set_value('rip_input', Main.HRU_Input([1.0,0,0,0], 0.01,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.2, np.zeros(4), 0, 0.0))
        model.set_value('Total_Elevationbands', 4)
        model.set_value('Elevation_Percentage', [0.16,0.46,0.33,0.05])
        model.set_value('bare_storage', Main.Storages(0,np.zeros(4),np.zeros(4),np.zeros(4),0))
        model.set_value('forest_storage', Main.Storages(0,np.zeros(4),np.zeros(4),np.zeros(4),0))
        model.set_value('grass_storage', Main.Storages(0,np.zeros(4),np.zeros(4),np.zeros(4),0))
        model.set_value('rip_storage', Main.Storages(0,np.zeros(4),np.zeros(4),np.zeros(4),0))


        # ## Running the model

        Discharge = []
        timestamp = []
        while (model.get_value_ptr('Current_Date') < (model.get_value_ptr('Date')[-1])):    
            model.update()
            timestamp.append(model.get_value_ptr('Current_Date'))
            Discharge.append(model.get_value_ptr('Discharge'))


        # ## Analysis


        simulated_discharge_df = pd.DataFrame(
            {f'simulation_{i+1}': Discharge},
            index=pd.to_datetime(timestamp)
        )
        # combined_discharge = pd.merge(simulated_discharge_df, forcing['streamflow'], left_index=True, right_index=True)
        # combined_discharge
        df = pd.merge(df, simulated_discharge_df, left_index=True, right_index=True)
        model.finalize()
    return df




def HBVmountain_simulation(ntimes):
    ob_list = []
    params_list = []
    forcing = pd.read_csv('Data/BigCreek/forcing_bigrockcreek.csv', index_col=[0], parse_dates=True)
    pd.to_datetime(forcing.index);
    forcing = forcing.reset_index(level=0)
    for i in range(len(forcing)):
        forcing['Date'][i] = forcing['Date'][i].date()
    forcing.set_index('Date', inplace=True)
    
    for i in range(ntimes):
        model = BMI_HBVmountain()

        config_file = model.setup()

        model.initialize(config_file)

        bare_parameters, forest_parameters, grass_parameters, rip_parameters,  slow_parameters, parameters_array =  Main.parameter_selection()
        model.set_value('bare_parameters', bare_parameters)
        model.set_value('forest_parameters', forest_parameters)
        model.set_value('grass_parameters', grass_parameters)
        model.set_value('rip_parameters', rip_parameters)
        model.set_value('slow_parameters', slow_parameters)

        model.set_value('Temperature', (forcing['monthly_temp_cru'].values).reshape(len(forcing),1))
        model.set_value('Precipitation', (forcing['prec_REGEN'].values).reshape(len(forcing),1))

        model.set_value('Date', list(forcing.index.values))
        model.set_value('Current_Date', forcing.index.values[0])

        model.set_value('Elevation', Main.Elevations(270, 1250, 2600, 1250, 1250))
        model.set_value('Glacier', [0.0, 0.0, 0.0, 0.0])
        model.set_value('Sunhours', [10.18, 10.90, 12.0, 13.10, 14.0, 14.45, 14.08, 13.31, 12.24, 11.25, 10.31, 9.85])
        model.set_value('bare_input', Main.HRU_Input([0,0,0,0], 0, np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.2, np.zeros(4), 0, 0.0))
        model.set_value('forest_input', Main.HRU_Input([0.0,0.20,0.73,0.07], 0.64,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.2, np.zeros(4), 0, 0.0))
        model.set_value('grass_input', Main.HRU_Input([0.1,0.8,0.1,0.0], 0.35,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.2, np.zeros(4), 0, 0.0))
        model.set_value('rip_input', Main.HRU_Input([1.0,0,0,0], 0.01,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.2, np.zeros(4), 0, 0.0))
        model.set_value('Total_Elevationbands', 4)
        model.set_value('Elevation_Percentage', [0.16,0.46,0.33,0.05])
        model.set_value('bare_storage', Main.Storages(0,np.zeros(4),np.zeros(4),np.zeros(4),0))
        model.set_value('forest_storage', Main.Storages(0,np.zeros(4),np.zeros(4),np.zeros(4),0))
        model.set_value('grass_storage', Main.Storages(0,np.zeros(4),np.zeros(4),np.zeros(4),0))
        model.set_value('rip_storage', Main.Storages(0,np.zeros(4),np.zeros(4),np.zeros(4),0))

        Discharge = []
        timestamp = []
        while (model.get_value_ptr('Current_Date') < (datetime.date(1996, 1, 1))):  
            model.update()
            timestamp.append(model.get_value_ptr('Current_Date'))
            Discharge.append(model.get_value_ptr('Discharge'))

        simulated_discharge_df = pd.DataFrame(
            {'simulation': Discharge},
            index=pd.to_datetime(timestamp)
        )
        combined_discharge = pd.merge(simulated_discharge_df, forcing['streamflow'], left_index=True, right_index=True)
        
        objective_function = multi_objective(combined_discharge['simulation'].values[1100:-1], combined_discharge['streamflow'].values[1100:-1], forcing['prec_REGEN'].values[1100:-1])
        ob_list.append(objective_function)
        params_list.append([bare_parameters, forest_parameters, grass_parameters, rip_parameters,  slow_parameters])
        
    bare_paramset = forest_paramset = grass_paramset = rip_paramset = pd.DataFrame(columns=['ED', 'NSE', 'logNSE', 'NSEfdc', NSErunoff, 'beta', 'Ce', 'Drainagecapacity', 'Interceptionstoragecapacity',
                              'Kf', 'Meltfactor', 'Mm', 'Ratio_Pref', 'Soilstoragecapacity', 'Temp_Thresh'])
    slow_paramset = pd.DataFrame(columns=['NSE', 'Ks', 'Ratio_Riparian'])

    for i in range(len(NSElist)):
        bare_paramset.loc[i] = [ob_list[i][0], ob_list[i][1], ob_list[i][2], ob_list[i][3], ob_list[i][4], params_list[i][0].beta, params_list[i][0].Ce, params_list[i][0].Drainagecapacity, 
                                params_list[i][0].Interceptionstoragecapacity,
                                     params_list[i][0].Kf, params_list[i][0].Meltfactor, params_list[i][0].Mm, params_list[i][0].Ratio_Pref, 
                                params_list[i][0].Soilstoragecapacity, params_list[i][0].Temp_Thresh]
        forest_paramset.loc[i] = [ob_list[i][0], ob_list[i][1], ob_list[i][2], ob_list[i][3], ob_list[i][4], params_list[i][1].beta, params_list[i][1].Ce, params_list[i][1].Drainagecapacity, 
                                  params_list[i][1].Interceptionstoragecapacity, params_list[i][1].Kf, params_list[i][1].Meltfactor,
                                  params_list[i][1].Mm, params_list[i][1].Ratio_Pref, params_list[i][1].Soilstoragecapacity, params_list[i][1].Temp_Thresh]
        grass_paramset.loc[i] = [ob_list[i][0], ob_list[i][1], ob_list[i][2], ob_list[i][3], ob_list[i][4], params_list[i][2].beta, params_list[i][2].Ce, params_list[i][2].Drainagecapacity,
                                 params_list[i][2].Interceptionstoragecapacity, params_list[i][2].Kf, params_list[i][2].Meltfactor, 
                                 params_list[i][0].Mm, params_list[i][2].Ratio_Pref, params_list[i][2].Soilstoragecapacity,
                                    params_list[i][2].Temp_Thresh]
        rip_paramset.loc[i] = [ob_list[i][0], ob_list[i][1], ob_list[i][2], ob_list[i][3], ob_list[i][4], params_list[i][3].beta, params_list[i][3].Ce, params_list[i][3].Drainagecapacity, 
                               params_list[i][3].Interceptionstoragecapacity,
        params_list[i][3].Kf, params_list[i][3].Meltfactor, params_list[i][3].Mm, 
                               params_list[i][3].Ratio_Pref, params_list[i][3].Soilstoragecapacity, params_list[i][3].Temp_Thresh]
        slow_paramset.loc[i] = [ob_list[i][0], ob_list[i][1], ob_list[i][2], ob_list[i][3], ob_list[i][4], params_list[i][4].Ks, params_list[i][4].Ratio_Riparian]
    return bare_paramset, forest_paramset, grass_paramset, rip_paramset, slow_paramset

def main():
    n_samples = int(sys.argv[1])
    bare_paramset, forest_paramset, grass_paramset, rip_paramset, slow_paramset = HBVmountain_simulation(n_samples)
    

    bare_name = 'bare_paramsets.csv'
    forest_name = 'forest_paramsets.csv'
    grass_name = 'grass_paramsets.csv'
    rip_name = 'rip_paramsets.csv'
    slow_name = 'slow_paramsets.csv'

    outdir = './output'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    bare_fullname = os.path.join(outdir, bare_name)
    forest_fullname = os.path.join(outdir, forest_name) 
    grass_fullname = os.path.join(outdir, grass_name) 
    rip_fullname = os.path.join(outdir, rip_name) 
    slow_fullname = os.path.join(outdir, slow_name) 


    bare_paramset.to_csv(bare_fullname)
    forest_paramset.to_csv(forest_fullname)
    grass_paramset.to_csv(grass_fullname)
    rip_paramset.to_csv(rip_fullname)
    slow_paramset.to_csv(slow_fullname)

if __name__ == '__main__':
    main()