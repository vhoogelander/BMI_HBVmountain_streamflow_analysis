"""Marrmot diagnostic."""
import logging
from pathlib import Path

import iris
import scipy.io as sio

from esmvalcore import preprocessor as preproc
from esmvaltool.diag_scripts.hydrology.derive_evspsblpot import debruin_pet
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic)

logger = logging.getLogger(Path(__file__).name)


def create_provenance_record():
    """Create a provenance record."""
    record = {
        'caption': "Forcings for the HBVmountain hydrological model.",
        'domains': ['global'],
        'authors': [
            'kalverla_peter',
            'camphuijsen_jaro',
            'alidoost_sarah',   #hoogelander_vincent
        ],
        'projects': [
            'ewatercycle',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': [],
    }
    return record


def get_input_cubes(metadata):
    """Return a dict with all (preprocessed) input files."""
    provenance = create_provenance_record()
    all_vars = {}
    for attributes in metadata:
        short_name = attributes['short_name']
        if short_name in all_vars:
            raise ValueError(
                f"Multiple input files found for variable '{short_name}'.")
        filename = attributes['filename']
        logger.info("Loading variable %s", short_name)
        all_vars[short_name] = iris.load_cube(filename)
        provenance['ancestors'].append(filename)
    return all_vars, provenance


def _get_extra_info(cube):
    """Get start/end time and origin of cube.

    Get the start and end times as an array with length 6
    and get latitude and longitude as an array with length 2
    """
    coord = cube.coord('time')
    time_start_end = []
    for index in 0, -1:
        time_val = coord.cell(index).point
        time_val = [
            time_val.year,
            time_val.month,
            time_val.day,
            time_val.hour,
            time_val.minute,
            time_val.second,
        ]
        time_val = [float(time) for time in time_val]
        time_start_end.append(time_val)

    # Add data_origin
    lat_lon = [
        cube.coord(name).points[0] for name in ('latitude', 'longitude')
    ]
    return time_start_end, lat_lon


def _shift_era5_time_coordinate(cube):
    """Shift instantaneous variables 30 minutes forward in time.

    After this shift, as an example:
    time format [1990, 1, 1, 11, 30, 0] will be [1990, 1, 1, 12, 0, 0].
    For aggregated variables, already time format is [1990, 1, 1, 12, 0, 0].
    """
    time = cube.coord(axis='T')
    time.points = time.points + 30 / (24 * 60)
    time.bounds = None
    time.guess_bounds()
    return cube

def save(cubes, dataset, exp, provenance, cfg):
    """Save cubes to file.

    Output format: "HBVmountain_local_forcing_ERA5_Meuse_1990_2018_scenario.nc"
    """
    time_coord = cubes[0].coord('time')
    start_year = time_coord.cell(0).point.year
    end_year = time_coord.cell(-1).point.year
    if dataset == "ERA5":
        basename = '_'.join([
            'HBVmountain',
            dataset,
            cfg['basin'],
            str(start_year),
            str(end_year),
        ])
    else:
        basename = '_'.join([
            'HBVmountain',
            dataset,
            cfg['basin'],
            str(start_year),
            str(end_year),
            exp,
        ])        
    output_file = get_diagnostic_filename(basename, cfg)
    logger.info("Saving cubes to file %s", output_file)
    iris.save(cubes, output_file, fill_value=-9999)

    # Store provenance
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(output_file, provenance)

def main(cfg):
    """Process data for use as input to the HBVmountain hydrological model.

    These variables are needed in all_vars:
    tas (air_temperature)
    pr (precipitation_flux)
    
    """
    input_metadata = cfg['input_data'].values()
    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        all_vars, provenance = get_input_cubes(metadata)
        if dataset != 'ERA5':
            exp = metadata[0]['exp']
        else:
            exp = ''
        # Fix time coordinate of ERA5 instantaneous variables
        if dataset == 'ERA5':

            _shift_era5_time_coordinate(all_vars['tas'])


        logger.info("Processing variable tas")
        temp = preproc.area_statistics(all_vars['tas'], operator='mean')


        logger.info("Processing variable pr")
        precip = preproc.area_statistics(all_vars['pr'], operator='mean')


        # Get the start and end times and latitude longitude
        time_start_end, lat_lon = _get_extra_info(temp)

        # Get the start and end times and latitude longitude
        time_start_end, lat_lon = _get_extra_info(temp) 

        # make data structure
        # delta_t_days could also be extracted from the cube
        output_data = {
            'forcing': {
                'precip': precip.data,
                'temp': temp.data,
#                'pet': pet.data,
                'delta_t_days': float(1),
                'time_unit': 'month',
            },
            'time_start': time_start_end[0],
            'time_end': time_start_end[1],
            'data_origin': lat_lon,
        }

       # Save to netcdf structure
        if dataset != 'ERA5':
            basename = '_'.join([
                'HBVmountain',
                dataset,
                cfg['basin'],
                str(int(output_data['time_start'][0])),
                str(int(output_data['time_end'][0])),
            ])
        else:
            basename = '_'.join([
                'HBVmountain',
                dataset,
                cfg['basin'],
                str(int(output_data['time_start'][0])),
                str(int(output_data['time_end'][0])),
                exp,
            ])            
        output_name = get_diagnostic_filename(basename, cfg, extension='nc')
        cubes = iris.cube.CubeList([temp, precip])
        save(cubes, dataset, exp, provenance, cfg)

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
