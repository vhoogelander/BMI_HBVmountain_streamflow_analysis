from julia.api import LibJulia
api = LibJulia.load()
api.sysimage = "sys.so"
api.init_julia()

from mpi4py import MPI
from BMI_HBVmountain_Python import *
from Calibration import *
import os
import sys

####### Forcing and observation data ##########################

def HBVmountain_simulation(nsmallest):
    ob_list = []
    params_list = parameter_conversion('output', nsmallest=nsmallest)
    forcing = pd.read_csv('Data/BigCreek/forcing_bigrockcreek.csv', index_col=[0], parse_dates=True)
    pd.to_datetime(forcing.index);
    forcing = forcing.reset_index(level=0)
    for i in range(len(forcing)):
        forcing['Date'][i] = forcing['Date'][i].date()
    forcing.set_index('Date', inplace=True)
    
    for i in range(nsmallest):
        model = BMI_HBVmountain()

        config_file = model.setup()

        model.initialize(config_file)

        model.set_value('bare_parameters', params_list[i][0])
        model.set_value('forest_parameters', params_list[i][1])
        model.set_value('grass_parameters', params_list[i][2])
        model.set_value('rip_parameters', params_list[i][3])
        model.set_value('slow_parameters', params_list[i][4])

        model.set_value('Temperature', (forcing['temp_era5'].values).reshape(len(forcing),1))
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
        while (model.get_value_ptr('Current_Date') < (model.get_value_ptr('Date')[-1])):   
            model.update()
            timestamp.append(model.get_value_ptr('Current_Date'))
            Discharge.append(model.get_value_ptr('Discharge'))

        simulated_discharge_df = pd.DataFrame(
            {'simulation': Discharge},
            index=pd.to_datetime(timestamp)
        )
        combined_discharge = pd.merge(simulated_discharge_df, forcing['streamflow'], left_index=True, right_index=True)
        
        
        objective_function = multi_objective(combined_discharge.simulation.loc[combined_discharge.index.year >= 1999], combined_discharge.streamflow.loc[combined_discharge.index.year >= 1999], forcing.prec_REGEN)
        
        
        ob_list.append(objective_function)
       
        model.finalize()
        
    bare_paramset = pd.DataFrame(columns=['ED', 'NSE', 'logNSE', 'NSEfdc', 'NSErunoff', 'beta', 'Ce', 'Drainagecapacity', 'Interceptionstoragecapacity',
                              'Kf', 'Meltfactor', 'Mm', 'Ratio_Pref', 'Soilstoragecapacity', 'Temp_Thresh'])
    forest_paramset = pd.DataFrame(columns=['ED', 'NSE', 'logNSE', 'NSEfdc', 'NSErunoff', 'beta', 'Ce', 'Drainagecapacity', 'Interceptionstoragecapacity',
                              'Kf', 'Meltfactor', 'Mm', 'Ratio_Pref', 'Soilstoragecapacity', 'Temp_Thresh'])
    grass_paramset = pd.DataFrame(columns=['ED', 'NSE', 'logNSE', 'NSEfdc', 'NSErunoff', 'beta', 'Ce', 'Drainagecapacity', 'Interceptionstoragecapacity',
                              'Kf', 'Meltfactor', 'Mm', 'Ratio_Pref', 'Soilstoragecapacity', 'Temp_Thresh'])
    rip_paramset = pd.DataFrame(columns=['ED', 'NSE', 'logNSE', 'NSEfdc', 'NSErunoff', 'beta', 'Ce', 'Drainagecapacity', 'Interceptionstoragecapacity',
                              'Kf', 'Meltfactor', 'Mm', 'Ratio_Pref', 'Soilstoragecapacity', 'Temp_Thresh'])
    slow_paramset = pd.DataFrame(columns=['ED', 'NSE', 'logNSE', 'NSEfdc', 'NSErunoff', 'Ks', 'Ratio_Riparian'])

    for i in range(len(ob_list)):
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
