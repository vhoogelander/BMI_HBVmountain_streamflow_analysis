from pandas import read_csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
from datetime import timedelta, date
import netCDF4 as nc

def generate_forcing_from_NETCDF(forcing_netcdf):
    path = path_to_forcing_netCDF
    ds = nc.Dataset(path)

    prec = (forcing_netcdf['pr'][:])*86400
    temp = (forcing_netcdf['tas'][:]) -273
    days = (forcing_netcdf['time'][:])

    start = date(1850,1,1)      # Netcdf file counts days from this date
  
    
    preclist = []
    tlist = []
    date_list=[]
    for i in range(len(temp)):
        preclist.append(prec[i])
        tlist.append(temp[i])
        delta = timedelta(days[i])
        date_list.append(start+delta)
    

    forcing = pd.DataFrame(index=date_list)
    forcing['temp'] = tlist
    forcing['prec'] = preclist
    
    
    forcing.loc[forcing['prec'] > 500, 'prec'] = 0  #remove outliers     
    return forcing
