from BMI_HBVmountain_Python import *
from Calibration import *
import os
import sys

####### Forcing and observation data ##########################
def run_validation(calibration_results):
    forcing = nc.Dataset('Data/BigCreek/HBVmountain_ERA5_BigRockCreek_1986_2015.nc') #Catchment dependent
    ob_list = []
    params_list = []
    sim_list = []
    
    df = pd.DataFrame(index=generate_forcing_from_NETCDF(forcing).prec.index)
    for i in range(len(calibration_results)):
        parameters = calibration_results.iloc[i, -20:]
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
        while (model.get_value_ptr('Current_Date') < (datetime.date(2005, 12, 31))):  
            model.update()
            timestamp.append(model.get_value_ptr('Current_Date'))
            Discharge.append(model.get_value_ptr('Discharge'))

        simulated_discharge =  pd.DataFrame(
                {'streamflow': Discharge},
                index=pd.to_datetime(timestamp)
            )
        simulated_discharge_df = pd.DataFrame(
            {f'simulation_{i+1}': Discharge},
            index=pd.to_datetime(timestamp)
        )
        df = pd.merge(df, simulated_discharge_df, left_index=True, right_index=True)
        
        model.finalize()
        
        # Validation 
        precipitation = generate_forcing_from_NETCDF(forcing).prec
        area = 59.3 #km2 #Catchment dependent
        observation = pd.read_csv('Data/BigCreek/Discharge_BigRockCreek.csv', index_col=0).streamflow / (area * 1e6) * 1000 *86400 #Catchment dependent
        precipitation.index, simulated_discharge.index, observation.index = pd.to_datetime(precipitation.index), pd.to_datetime(simulated_discharge.index), pd.to_datetime(observation.index)
        mask = (observation.index >= simulated_discharge.index[0]) & (observation.index <= simulated_discharge.index[-1])
        observation = observation.loc[mask]
        
        
        objective_function = multi_objective(simulated_discharge.loc[simulated_discharge.index.year >= 1999].streamflow, observation.loc[observation.index.year >= 1999], precipitation)
        
        
        ob_list.append(objective_function)
        params_list.append(parameters)    

        paramset = pd.DataFrame(columns=['ED', 'NSE', 'logNSE', 'NSEfdc', 'NSErunoff', 'Temp_Thresh', 'Meltfactor', 'Mm', 'Ratio_Pref', 'Kf', 'Ks', 'Ce', 'Soilstoragecapacity_Bare', 'beta_Bare', 'Interceptioncapacity_Forest', 'Soilstoragecapacity_Forest', 'beta_Forest', 'Interceptioncapacity_Grass', 'Soilstoragecapacity_Grass', 'beta_Grass', 'Interceptioncapacity_Rip', 'Soilstoragecapacity_Rip', 'beta_Rip', 'Kf_Rip', 'Ratio_Riparian'])
        
        
    for i in range(len(ob_list)):
        paramset.loc[i] = [ob_list[i][0], ob_list[i][1], ob_list[i][2], ob_list[i][3], ob_list[i][4], params_list[i][0], params_list[i][1], params_list[i][2], params_list[i][3], params_list[i][4], params_list[i][5], params_list[i][6], params_list[i][7], params_list[i][8], params_list[i][9], params_list[i][10], params_list[i][11], params_list[i][12], params_list[i][13], params_list[i][14], params_list[i][15], params_list[i][16], params_list[i][17], params_list[i][18], params_list[i][19]]    
        

    return paramset, df















def main():
    comm = MPI.COMM_WORLD
    cpus = comm.Get_size()
    rank = comm.Get_rank()
    n_samples = int(sys.argv[1])

    if rank == 0:
        start_time = datetime.datetime.now()
        partitions = [ int(n_samples / cpus) ] * cpus
        counts = [ int(0) ] * cpus
    else:
        partitions = None
        counts = None

    partition_item = comm.scatter(partitions, root=0)
    bare_paramset_item, forest_paramset_item, grass_paramset_item, rip_paramset_item, slow_paramset_item, = HBVmountain_simulation(partition_item)
    bare_paramset = comm.gather(bare_paramset_item, root=0)
    forest_paramset = comm.gather(forest_paramset_item, root=0)
    grass_paramset = comm.gather(grass_paramset_item, root=0)
    rip_paramset = comm.gather(rip_paramset_item, root=0)
    slow_paramset = comm.gather(slow_paramset_item, root=0)
    #bare_paramset, forest_paramset, grass_paramset, rip_paramset, slow_paramset = HBVmountain_simulation(n_samples)
    if rank == 0:
        bare_paramset = pd.concat(bare_paramset)
        forest_paramset = pd.concat(forest_paramset)
        grass_paramset = pd.concat(grass_paramset)
        rip_paramset = pd.concat(rip_paramset)
        slow_paramset = pd.concat(slow_paramset)


        bare_name = 'bare_paramsets.csv'
        forest_name = 'forest_paramsets.csv'
        grass_name = 'grass_paramsets.csv'
        rip_name = 'rip_paramsets.csv'
        slow_name = 'slow_paramsets.csv'

        outdir = './output_validation'
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
