
from BMI_HBVmountain_Python import *
from Calibration_BRC import *
import netCDF4 as nc

forcing = pd.read_csv('Data/BigCreek/forcing_bigrockcreek.csv', index_col=[0], parse_dates=True)
pd.to_datetime(forcing.index);
forcing = forcing.reset_index(level=0)
for i in range(len(forcing)):
    forcing['Date'][i] = forcing['Date'][i].date()
forcing.set_index('Date', inplace=True)
forcing


# ## Select parameter sets

calibration = pd.read_csv('output/bare_paramsets.csv').sort_values(by=['ED']).iloc[: , :6]
calibration = calibration.loc[calibration['ED'].isin(calibration.ED.nsmallest(60))]
calibration.reset_index(inplace=True)

validation  = pd.read_csv('output_validation/forest_paramsets.csv').iloc[: , :6]
validation.rename(columns={"ED": "ED_val", "NSE": "NSE_val", "logNSE": "logNSE_val", "NSEfdc": "NSEfdc_val", "NSErunoff": "NSErunoff_val"}, inplace=True)

calval_paramsets =  pd.concat([calibration, validation], axis=1)

calval_paramsets.loc[calval_paramsets.ED_val <=0.4]


# ## Historical climate data

path = 'Data/BigCreek/HBVmountain_GFDL-CM4_BigCreek_1975_2005.nc'
ds = nc.Dataset(path)

prec = (ds['pr'][:])*86400
temp = (ds['tas'][:]) -273

preclist = []
tlist = []

for i in range(len(temp)):
    preclist.append(prec[i])
    tlist.append(temp[i])

base = datetime.datetime.strptime('1975-01-01', '%Y-%m-%d')
date_list = [(base + datetime.timedelta(days=x)).date() for x in range(int(ds['time'][:][-1] - ds['time'][:][0] + 1))]

hist_forcing = pd.DataFrame({'prec':preclist, 'temp':tlist}, index=date_list)

start_date = pd.to_datetime('1975-01-01').date()
end_date = pd.to_datetime('2004-12-31').date()
# Select DataFrame rows between two dates
mask = (hist_forcing.index >= start_date) & (hist_forcing.index <= end_date)
hist_forcing = hist_forcing.loc[mask]


# ## Future climate forcing

path = 'Data/BigCreek/HBVmountain_GFDL-CM4_BigCreek_2050_2100.nc'
ds = nc.Dataset(path)

prec = (ds['pr'][:])*86400
temp = (ds['tas'][:]) -273

preclist = []
tlist = []

for i in range(len(temp)):
    preclist.append(prec[i])
    tlist.append(temp[i])

base = datetime.datetime.strptime('2050-01-01', '%Y-%m-%d')
date_list = [(base + datetime.timedelta(days=x)).date() for x in range(int(ds['time'][:][-1] - ds['time'][:][0] + 1))]

future_forcing = pd.DataFrame({'prec':preclist, 'temp':tlist}, index=date_list)
start_date = pd.to_datetime('2070-01-01').date()
end_date = pd.to_datetime('2099-12-31').date()
# Select DataFrame rows between two dates
mask = (future_forcing.index >= start_date) & (future_forcing.index <= end_date)
future_forcing = future_forcing.loc[mask]

param_list = [parameter_conversion('output', nsmallest=60)[i] for i in calval_paramsets.loc[calval_paramsets.ED_val <=0.4].index]


def climate_simulations(prec, temp, param_list):
    # ## Setting up the model
    forcing = pd.DataFrame({'prec':prec, 'temp':temp}, index=prec.index)
    df = pd.DataFrame()
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
            {'simulation': f'simulation {i}'},
            index=pd.to_datetime(timestamp)
        )
        # combined_discharge = pd.merge(simulated_discharge_df, forcing['streamflow'], left_index=True, right_index=True)
        # combined_discharge
        df.append(simulated_discharge_df)
        model.finalize()
    return df