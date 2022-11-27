import julia
from julia import Main
from pandas import read_csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
from datetime import timedelta, date
import netCDF4 as nc
from rasterio.features import shapes
from shapely.geometry import shape
import geopandas as gp
from rasterio import features
import rasterio
import json
from rasterio.mask import mask
import numpy.ma as ma
from suntime import Sun, SunTimeException
import cma
import geopandas as gpd
import fiona
import time
from bmipy import Bmi
from warnings import filterwarnings
import os
from scipy import stats
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

from ruamel.yaml import YAML
def generate_forcing_from_NETCDF(forcing_netcdf):


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
    forcing.index = pd.to_datetime(forcing.index)
    
    forcing.loc[forcing['prec'] > 500, 'prec'] = 0  #remove outliers 
    freq = pd.offsets.Hour(5)
    if freq.is_year_start((forcing.index[0])) == False:
        start = forcing.index[0] + pd.offsets.YearBegin()
        forcing = forcing.loc[start:]
    if freq.is_year_end((forcing.index[-1])) == False:
        end = forcing.index[-1] - pd.offsets.YearBegin()
        forcing = forcing.loc[:end][0:-1]
    forcing.index = pd.to_datetime(forcing.index).date
    return forcing


def generate_array_from_raster(str_path_to_rasterfile, str_path_to_shapefile=None):
    if str_path_to_shapefile != None:
        with fiona.open(str_path_to_shapefile, "r") as shapefile:
            shapefile = [feature["geometry"] for feature in shapefile]

        with rasterio.open(str_path_to_rasterfile) as src:
            out_image, out_transform = rasterio.mask.mask(src, shapefile, crop=True)
#             out_image = out_image[out_image <97]
            out_image = ma.masked_values(out_image, 0)
            out_image = out_image.astype(np.float32)
            out_image = out_image.filled(np.nan)
            arr = np.expand_dims(out_image.flatten(), 0).T    
    else:
        raster = rasterio.open(str_path_to_rasterfile)
        arr = raster.read(masked=True)
        arr = np.expand_dims(arr.flatten(), 0).T    
    return arr

def get_elevations_from_raster(raster_array_flattened, nbands):
   
    n, bins, patches = plt.hist(raster_array_flattened, bins=nbands)
    plt.close()
    tot_pixels = n.sum()
#     raster_array_flattened = raster_array_flattened.filled(np.nan)
    av_elevation = round(np.nanmean(raster_array_flattened),9)
#     av_elevation = raster_array_flattened.mean()
    elevation_list = []
    for i in range(nbands):
        elevation_percentage = n[i] / tot_pixels
        elevation_list.append(elevation_percentage)
    min_elevation = np.min(bins)
    max_elevation = np.max(bins)
    return elevation_list, bins, min_elevation, max_elevation, av_elevation

def get_landuse_from_raster(raster_array_flattened):
    n, bins, patches = plt.hist(raster_array_flattened, bins=np.arange(96))
    plt.close()
    
    tot_pixels = n.sum()
    
    bare = np.sum(n[[12, 22, 22, 24, 31]]) / tot_pixels 
    forest = np.sum(n[40:46]) / tot_pixels 
    grass = (np.sum(n[50:83])+ n[21]) /tot_pixels 
    rip = np.sum(n[[11, 90, 94]]) /tot_pixels  
    
    

    landuse_list = [bare, forest, grass, rip] 
    return landuse_list

def generate_shapefiles_per_landuse(landuse, str_path_to_rasterfile, str_path_to_shapefile=None):  
        
    with rasterio.Env():
        if str_path_to_shapefile != None:
            with fiona.open(str_path_to_shapefile, crs='epsg:4326') as shapefile:
                shapefile = [feature["geometry"] for feature in shapefile]
            with rasterio.open(str_path_to_rasterfile) as src:
                image, out_transform = rasterio.mask.mask(src, shapefile, crop=True)
                image = ma.masked_values(image, 0)
                image = image.astype(np.float32)
                image = image.filled(np.nan)

        else:
            with rasterio.open(str_path_to_rasterfile) as src:
                image = src.read() # first band
        if landuse == 'glacier':
            image[(image == 12)] = -9999
        if landuse == 'bare':
            image[(image == 12) | (image == 22) | (image == 23)  | (image == 24) | (image == 31)] = -9999
        if landuse == 'forest':
            image[(image == 41) | (image == 42) | (image == 43)] = -9999
        if landuse == 'grass':
            image[(image == 21) | (image == 51) | (image == 52) | (image == 71) | (image == 72) | (image == 73) | (image == 74) | (image == 81) | (image == 82)] = -9999
        if landuse == 'rip':
            image[(image == 11) | (image == 90) | (image == 95)] = -9999


        if str_path_to_shapefile != None:
            results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) 
            in enumerate(
                shapes(image, mask=(image ==-9999), transform= out_transform)))
        else:
            results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) 
            in enumerate(
                shapes(image, mask=(image ==-9999), transform= src.transform)))            
    geoms = list(results)
    gpd_polygonized_raster  = gp.GeoDataFrame.from_features(geoms, crs='epsg:4326')
    return gpd_polygonized_raster

def clip_elevation_per_landuse(gpd_polygonized_raster, str_path_to_elevation, str_path_to_shapefile=None, export_tif=False):
    gdf = gpd_polygonized_raster.dissolve()
    coords = [json.loads(gdf.to_json())['features'][0]['geometry']] 
    if str_path_to_shapefile != None:

        with fiona.open(str_path_to_shapefile) as shapefile:
            shapefile = [feature["geometry"] for feature in shapefile]
        with rasterio.open(str_path_to_elevation) as src:
#             out_img, out_transform = rasterio.mask.mask(src, shapefile, crop=True)
            out_img, out_transform = rasterio.mask.mask(src, coords, crop=True)
            out_img = ma.masked_values(out_img, 0)
            out_img = out_img.astype(np.float32)
            out_img = out_img.filled(np.nan) 
    
    else:
        data = rasterio.open(str_path_to_elevation, shapefile=None)
        out_img, out_transform = mask(data, shapes=coords)
        out_img = np.expand_dims(out_img.flatten(), 0).T
        out_img = ma.masked_values(out_img, data.nodata)
    
    if export_tif == True:
        out_meta = data.meta.copy()
        with rasterio.open('elevationlanduse.tif', "w", **out_meta) as dest:
            dest.write(out_img)
    
    return out_img 

def generate_landuse_per_elevation(str_path_to_elevation, str_path_to_landuse, str_path_to_shapefile=None, nbands=4):
    elevation_array = generate_array_from_raster(str_path_to_elevation, str_path_to_shapefile)
    landuse_array = generate_array_from_raster(str_path_to_landuse, str_path_to_shapefile)
    
    elevation_list, bins, min_elevation, max_elevation, av_elevation = get_elevations_from_raster(elevation_array, nbands=nbands)
    landuse_catchment = get_landuse_from_raster(landuse_array)
    
    tot_elevations = [[elevation_list, landuse_catchment, (max_elevation-min_elevation)/nbands, min_elevation, max_elevation, av_elevation]]
    landuse = ['glacier', 'bare', 'forest', 'grass', 'rip']
    for x in landuse:
        gpd_polygonized_raster = generate_shapefiles_per_landuse(x, str_path_to_landuse, str_path_to_shapefile)
        if not gpd_polygonized_raster.empty:
            arr = clip_elevation_per_landuse(gpd_polygonized_raster, str_path_to_elevation, str_path_to_shapefile)
            raster_array_flattened = np.expand_dims(arr.flatten(), 0).T
            elevation_list, bins, min_elevation, max_elevation, av_elevation = get_elevations_from_raster(raster_array_flattened, nbands=nbands)
        if gpd_polygonized_raster.empty:
            elevation_list = list(np.zeros(4))
        if x == 'glacier':
            tot_elevations.append([elevation_list, list(np.around(np.array(elevation_list)*landuse_catchment[0], 9))]) 
        else:
            tot_elevations.append(elevation_list)

#         tot_elevations.append([elevation_list])
    return tot_elevations

def linear_scaling_correction(str_path_calibration_forcing, str_path_input_hist_forcing, str_path_input_forcing):
    calibration_forcing = nc.Dataset(str_path_calibration_forcing)
    input_hist_forcing = nc.Dataset(str_path_input_hist_forcing)
    input_forcing = nc.Dataset(str_path_input_forcing)

    calibration_forcing  = generate_forcing_from_NETCDF(calibration_forcing)
    calibration_forcing.index = pd.to_datetime(calibration_forcing.index)
    calibration_forcing = calibration_forcing.loc[calibration_forcing.index.year <= calibration_forcing.index.year[-1]]

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
    
def get_monthly_sunhours(lat, lon):
    sunhours = []
    sun = Sun(lat, lon)
    for i in range(1,13):
        timedelta = sun.get_local_sunset_time(datetime.date(2005, i, 14)) - sun.get_local_sunrise_time(datetime.date(2005, i, 14))
        sunhours.append(np.round(timedelta.seconds / 3600,2))
    return sunhours

def select_cmip6_models(str_path_to_models, model_list=None, ignore_models=None, select_all=False):
    model_file = yaml.safe_load(open(str_path_to_models))
            
    if ignore_models != None:
        selected_CMIP6_models = dict((k, model_file[k]) for k in model_file
           if k not in ignore_models)
        selected_CMIP6_models = tuple(selected_CMIP6_models.values()) 
    if select_all == True and model_list == None:
        selected_CMIP6_models = dict((k, model_file[k]) for k in model_file.keys()
               if k in model_file)
        selected_CMIP6_models = tuple(selected_CMIP6_models.values()) 
    if model_list != None:
        selected_CMIP6_models = dict((k, model_file[k]) for k in model_list
               if k in model_file)
        selected_CMIP6_models = tuple(selected_CMIP6_models.values())
    return selected_CMIP6_models

def run_recipe_HBVmountain(recipe, catchment_name, scenario, CMIP6_models, start_year, end_year, catchment_id=None):
    recipe.data["preprocessors"]["daily"]["extract_shape"]["shapefile"] = os.path.join(catchment_name, catchment_name + '.' + 'shp') #shape
    recipe.data["diagnostics"]["diagnostic_daily"]["scripts"]["script"][
        "basin"
    ] = catchment_name #basin
    
    recipe.data["diagnostics"]["diagnostic_daily"]["additional_datasets"] = [
    *CMIP6_models
]
    variables = recipe.data["diagnostics"]["diagnostic_daily"]["variables"]
    var_names = "tas", "pr"
    
    for var_name in var_names:
        variables[var_name]["exp"] = scenario
        
        
    startyear = start_year #get_time(start_time).year
    for var_name in var_names:
        variables[var_name]["start_year"] = startyear

    endyear = end_year # get_time(end_time).year
    for var_name in var_names:
        variables[var_name]["end_year"] = endyear
    recipe_output = recipe.run()
    return recipe_output