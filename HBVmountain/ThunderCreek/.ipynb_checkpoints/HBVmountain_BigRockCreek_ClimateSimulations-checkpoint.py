#!/usr/bin/env python
# coding: utf-8

# In[1]:


from BMI_HBVmountain_Python import *
from Calibration_BRC import *
import netCDF4 as nc


# ## Load forcing and observation data

# In[10]:


forcing = pd.read_csv('Data/BigCreek/forcing_bigrockcreek.csv', index_col=[0], parse_dates=True)
pd.to_datetime(forcing.index);
forcing = forcing.reset_index(level=0)
for i in range(len(forcing)):
    forcing['Date'][i] = forcing['Date'][i].date()
forcing.set_index('Date', inplace=True)
forcing


# ## Select parameter sets

# In[282]:


calibration = pd.read_csv('output/bare_paramsets.csv').sort_values(by=['ED']).iloc[: , :6]
calibration = calibration.loc[calibration['ED'].isin(calibration.ED.nsmallest(60))]
calibration.reset_index(inplace=True)

validation  = pd.read_csv('output_validation/forest_paramsets.csv').iloc[: , :6]
validation.rename(columns={"ED": "ED_val", "NSE": "NSE_val", "logNSE": "logNSE_val", "NSEfdc": "NSEfdc_val", "NSErunoff": "NSErunoff_val"}, inplace=True)

calval_paramsets =  pd.concat([calibration, validation], axis=1)

calval_paramsets.loc[calval_paramsets.ED_val <=0.4]


# ## Historical climate data

# In[340]:


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

# In[258]:


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


# In[259]:


start_date = pd.to_datetime('2070-01-01').date()
end_date = pd.to_datetime('2099-12-31').date()
# Select DataFrame rows between two dates
mask = (future_forcing.index >= start_date) & (future_forcing.index <= end_date)
future_forcing = future_forcing.loc[mask]

future_forcing


# ## Setting up the model

# In[387]:


model = BMI_HBVmountain()


# In[388]:


config_file = model.setup()


# In[389]:


model.initialize(config_file)


# ### Parameters

# In[390]:


# param_list = parameter_conversion('output', nsmallest=60).iloc[calval_paramsets.loc[calval_paramsets.ED_val <=0.4]]
param_list = [parameter_conversion('output', nsmallest=60)[i] for i in calval_paramsets.loc[calval_paramsets.ED_val <=0.4].index]
model.set_value('bare_parameters', param_list[0][0])
model.set_value('forest_parameters', param_list[0][1])
model.set_value('grass_parameters', param_list[0][2])
model.set_value('rip_parameters', param_list[0][3])
model.set_value('slow_parameters', param_list[0][4])


# ## Historical forcing

# In[391]:


model.set_value('Temperature', (hist_forcing['temp'].values).reshape(len(hist_forcing),1))
model.set_value('Precipitation', (hist_forcing['prec'].values).reshape(len(hist_forcing),1))

model.set_value('Date', list(hist_forcing.index.values))
model.set_value('Current_Date', hist_forcing.index.values[0])


# ### Future forcings

# In[380]:


model.set_value('Temperature', (future_forcing['temp'].values).reshape(len(future_forcing),1))
model.set_value('Precipitation', (future_forcing['prec'].values).reshape(len(future_forcing),1))

model.set_value('Date', list(future_forcing.index.values))
model.set_value('Current_Date', future_forcing.index.values[0])


# ### Initial settings

# In[392]:


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

# In[393]:


Discharge = []
timestamp = []
# while (model.get_value_ptr('Current_Date') < (model.get_value_ptr('Date')[-1])):  
while (model.get_value_ptr('Current_Date') < (datetime.date(2004, 12, 31))):  
    model.update()
    timestamp.append(model.get_value_ptr('Current_Date'))
    Discharge.append(model.get_value_ptr('Discharge'))


# ## Analysis

# In[394]:


simulated_discharge_df = pd.DataFrame(
    {'simulation': Discharge},
    index=pd.to_datetime(timestamp)
)
# combined_discharge = pd.merge(simulated_discharge_df, forcing['streamflow'], left_index=True, right_index=True)
# combined_discharge


# In[395]:


plt.figure(figsize=[16,8])
simulated_discharge_df.simulation.plot(label='Simulation')
plt.legend();

model.get_value_ptr('Current_Date')


# In[396]:


rolling = simulated_discharge_df.simulation.rolling(window=7).mean()
#     df['MA'] = df.rolling(window=5).mean()
histdoy = []
histmagn = []
for i in range(1975, 2004):
#     rolling_year = pd.DataFrame(rolling.loc[simulated_discharge_df.index.year == i])
    rolling_year = pd.DataFrame(rolling.loc[(simulated_discharge_df.index >= datetime.datetime(i, 3, 1)) & (simulated_discharge_df.index < datetime.datetime(i, 11, 1))])
    rolling_year.reset_index(inplace=True)
    rolling_year.rename(columns={"index":"Date"}, inplace=True)
    rolling_year['cs'] = rolling_year.simulation.cumsum()
#     plt.figure()
#     plt.plot(rolling_year.Date, rolling_year.cs)
    plt.title(f'{rolling_year.loc[rolling_year.cs >= (0.75*rolling_year.cs.max())].index[0]}, min_flow= {rolling_year.cs.min():.3f}')
    histdoy.append(rolling_year.loc[rolling_year.cs >= (0.75*rolling_year.cs.max())].index[0])
    histmagn.append(rolling_year.cs.min())
#     doy.append(rolling_year.simulation.idxmin())
#     plt.figure()
#     plt.plot(rolling_year.Date, rolling_year.simulation)
#     plt.title(rolling_year.simulation.idxmin())      
hist = [121, 122, 123, 123, 122, 121, 118, 127, 125, 130, 118, 117, 124, 125]
histmagn = [0.7599267848259975, 0.2813994807338532,1.0771969104804655,4.009510870996008,0.6676913212599985,1.3357459933440057,1.0473084360001281,0.6438389960903111,4.008755979821798,0.18316181441137466,0.6558363343211738,0.6020336764498682,0.11030228747160672,0.676352208858607]
np.mean(doy), np.mean(magn), np.mean(histdoy), np.mean(histmagn)


# In[242]:


np.mean(magn)


# In[ ]:


rolling_year


# ## Clean up

# In[ ]:


model.finalize()


# In[201]:


len(model.get_value_ptr('Temperature')


# In[149]:


model.get_value_ptr('Temperature')[0:20]


# In[64]:


directory = 'output'

bare_parametersets = pd.read_csv(os.path.join(directory, 'bare_paramsets.csv')).sort_values(by=['ED'])
# bare_parametersets = bare_parametersets.loc[bare_parametersets['ED'].isin(bare_parametersets.ED.nsmallest(25))]
forest_parametersets = pd.read_csv(os.path.join(directory, 'forest_paramsets.csv')).loc[bare_parametersets.index].sort_values(by=['ED'])
grass_parametersets = pd.read_csv(os.path.join(directory, 'grass_paramsets.csv')).loc[bare_parametersets.index].sort_values(by=['ED'])
rip_parametersets = pd.read_csv(os.path.join(directory, 'rip_paramsets.csv')).loc[bare_parametersets.index].sort_values(by=['ED'])
slow_parametersets = pd.read_csv(os.path.join(directory, 'slow_paramsets.csv')).loc[bare_parametersets.index].sort_values(by=['ED'])


# In[65]:


for i in np.arange(6, 16, 1):
    plt.figure()
    plt.xlabel='ED'
    for j in range(len(forest_parametersets)):
        plt.plot(forest_parametersets.iloc[j, 1], forest_parametersets.iloc[j, i], '.', color='blue')


# In[56]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




