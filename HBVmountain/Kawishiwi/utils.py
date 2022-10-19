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
    
    
    forcing.loc[forcing['prec'] > 500, 'prec'] = 0  #remove outliers     
    return forcing

def generate_array_from_raster(str_path_to_rasterfile):
    raster = rasterio.open(str_path_to_rasterfile)
    arr = raster.read(masked=True)
    arr = np.expand_dims(arr.flatten(), 0).T
    
    return arr

def get_elevations_from_raster(raster_array_flattened, nbands=4):
   
    n, bins, patches = plt.hist(raster_array_flattened, bins=nbands)
    plt.close()
    tot_pixels = n.sum()
    
    elevation_list = []
    for i in range(nbands):
        elevation_percentage = n[i] / tot_pixels
        elevation_list.append(round(elevation_percentage, 3))
    min_elevation = np.min(bins)
    max_elevation = np.max(bins)
    return elevation_list, bins, min_elevation, max_elevation

def get_landuse_from_raster(raster_array_flattened):
    n, bins, patches = plt.hist(raster_array_flattened, bins=np.arange(96))
    plt.close()
    
    tot_pixels = n.sum()
    
    bare = np.sum(n[[11, 21, 22, 23, 30]]) / tot_pixels 
    forest = np.sum(n[40:44]) / tot_pixels 
    grass = (np.sum(n[50:82])+ n[20]) /tot_pixels 
    rip = np.sum(n[[10, 89, 94]]) /tot_pixels  
    
    landuse_list = [bare, forest, grass, rip] 
    return landuse_list

def generate_shapefiles_per_landuse(landuse, str_path_to_nlcd):
    with rasterio.Env():
        with rasterio.open(str_path_to_nlcd) as src:
            image = src.read(masked=True) # first band
                    
            if landuse == 'bare':
                image[(image == 12) | (image == 22) | (image == 23)  | (image == 24) | (image == 31)] = -9999
            if landuse == 'forest':
                image[(image == 41) | (image == 42) | (image == 43)] = -9999
            if landuse == 'grass':
                image[(image == 21) | (image == 51) | (image == 52) | (image == 71) | (image == 72) | (image == 73) | (image == 74) | (image == 81) | (image == 82)] = -9999
            if landuse == 'rip':
                image[(image == 11) | (image == 90) | (image == 95)] = -9999
                
            results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) 
            in enumerate(
                shapes(image, mask=(image ==-9999), transform=src.transform)))
    geoms = list(results)
    gpd_polygonized_raster  = gp.GeoDataFrame.from_features(geoms, crs='epsg:4326')
    return gpd_polygonized_raster

def clip_elevation_per_landuse(gpd_polygonized_raster, str_path_to_elevation, export_tif=False):
    gdf = gpd_polygonized_raster.dissolve()
    coords = [json.loads(gdf.to_json())['features'][0]['geometry']] #parse features from GeoDataFrame in such a manner that rasterio wants them


    data = rasterio.open(str_path_to_elevation, masked=True)
    out_img, out_transform = mask(data, shapes=coords)
    out_img = np.expand_dims(out_img.flatten(), 0).T
    out_img = ma.masked_values(out_img, data.nodata)
    if export_tif == True:
        out_meta = data.meta.copy()
        with rasterio.open('elevationlanduse.tif', "w", **out_meta) as dest:
            dest.write(out_img)
    
    return out_img 

def generate_landuse_per_elevation(str_path_to_elevation, str_path_to_landuse):
    elevation_array = generate_array_from_raster(str_path_to_elevation)
    landuse_array = generate_array_from_raster(str_path_to_landuse)
    
    elevation_list, bins, min_elevation, max_elevation = get_elevations_from_raster(elevation_array, nbands=4)
    landuse_catchment = get_landuse_from_raster(landuse_array)
    
    tot_elevations = [elevation_list]
    landuse = ['bare', 'forest', 'grass', 'rip']
    for x in landuse:
        gpd_polygonized_raster = generate_shapefiles_per_landuse(x, str_path_to_landuse)
        arr = clip_elevation_per_landuse(gpd_polygonized_raster, str_path_to_elevation)
        raster_array_flattened = np.expand_dims(arr.flatten(), 0).T
        elevation_list, bins, min_elevation, max_elevation = get_elevations_from_raster(raster_array_flattened, nbands=4)
        tot_elevations.append([elevation_list, min_elevation])
    return tot_elevations