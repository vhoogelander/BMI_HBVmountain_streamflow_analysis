#!/usr/bin/env python
# coding: utf-8

# In[1]:


from BMI_HBVmountain_Python import *
from Calibration_BRC import *
import os
import netCDF4 as nc

import cma
import hashlib
import time
import warnings
from pathlib import Path
# from cartopy.io import shapereader
from cftime import date2num
from matplotlib.ticker import MaxNLocator
import logging
logging.basicConfig(level=logging.WARNING)
warnings.filterwarnings("ignore", category=UserWarning)
from cma.fitness_transformations import EvalParallel2

# import fiona
# from cartopy.io import shapereader
# from ruamel.yaml import YAML
# import ewatercycle.models
# from ewatercycle.config import CFG
# from ewatercycle.forcing import load_foreign
# from ewatercycle.observation.grdc import get_grdc_data
# from ewatercycle.util import get_time


# In[2]:


POPSIZE = 2  # it can be equal to number of available cores * 0.75
MAXITER = 2 


#  ### Load forcing data

# In[3]:


forcing = pd.read_csv('Data/ThunderCreek/forcing_thundercreek.csv', index_col=[0], parse_dates=True)
pd.to_datetime(forcing.index);
forcing = forcing.reset_index(level=0)
for i in range(len(forcing)):
    forcing['time'][i] = forcing['time'][i].date()
forcing.set_index('time', inplace=True)


# ### Functions

# In[4]:


def run_model_cma(parameters):
    """Setup and run model."""
    model = BMI_HBVmountain()
    

    # Generate a random tag for work_dir name
    # Each work_dir is used by one worker in multiprocessing
#    tag = hashlib.md5(str(parameters).encode('utf-8')).hexdigest()
#    work_dir = CFG['output_dir'] / f"marrmot-m14_{tag}"
    
    
    if len(parameters) ==5:
        cfg_file = model.setup(bare_parameters=  Parameters(parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3], parameters[0][4],
                                                                     parameters[0][5], parameters[0][6], parameters[0][7], parameters[0][8], parameters[0][9]),
                                        forest_parameters=Parameters(parameters[1][0], parameters[1][1], parameters[1][2], parameters[1][3], parameters[1][4],
                                                                     parameters[1][5], parameters[1][6], parameters[1][7], parameters[1][8], parameters[1][9]), 
                                        grass_parameters= Parameters(parameters[2][0], parameters[2][1], parameters[2][2], parameters[2][3], parameters[2][4],
                                                                     parameters[2][5], parameters[2][6], parameters[2][7], parameters[2][8], parameters[2][9]),
                                        rip_parameters=   Parameters(parameters[3][0], parameters[3][1], parameters[3][2], parameters[3][3], parameters[3][4],
                                                                     parameters[3][5], parameters[3][6], parameters[3][7], parameters[3][8], parameters[3][9]), 
                                        slow_parameters=  Slow_Paramters(parameters[4][0], parameters[4][1])
        )
    else:
        cfg_file = model.setup(bare_parameters=  Parameters(parameters[8], parameters[6], 0, parameters[0], parameters[4],
                                                                     parameters[1], parameters[2], parameters[3], parameters[7], parameters[0]),
                                        forest_parameters=Parameters(parameters[11], parameters[6], 0, parameters[9], parameters[4],
                                                                     parameters[1], parameters[2], parameters[3], parameters[10], parameters[0]),
                                        grass_parameters= Parameters(parameters[14], parameters[6], 0, parameters[12], parameters[4],
                                                                     parameters[1], parameters[2], parameters[3], parameters[13], parameters[0]),
                                        rip_parameters=   Parameters(parameters[17], parameters[6], 0, parameters[15], parameters[18],
                                                                     parameters[1], parameters[2], parameters[3], parameters[16], parameters[0]),
                                        slow_parameters=  Slow_Paramters(parameters[5], parameters[19]))        
    
    # Initialize the model inside the docker container through GRPC4bmi
    model.initialize(cfg_file)
    
    model.set_value('Temperature', (forcing['temp_era5'].values).reshape(len(forcing),1))
    model.set_value('Precipitation', (forcing['prec_era5'].values).reshape(len(forcing),1))
    model.set_value('Date', list(forcing.index.values))
    model.set_value('Current_Date', forcing.index.values[0])  
    
    model.set_value('Elevation', Elevations(500, 500, 2500, 1500, 1500)) # 4 Elevation bands
    model.set_value('Glacier', [0.0, 0.0, 0.0, 0.6]) # 60% of the area glacier at highest elevation band
    model.set_value('Sunhours', [8.87, 10.30, 11.88, 13.65, 15.13, 15.97, 15.58, 14.25, 12.62, 11.87, 9.28, 8.45]) #Seattle
    model.set_value('bare_input', HRU_Input([0.0,0.0,0.3,0.7], 0.32, [0.0, 0.0, 0.0, 0.6], [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
    model.set_value('forest_input', HRU_Input([0.0,0.7,0.3,0.0], 0.45,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
    model.set_value('grass_input', HRU_Input([0.7,0.3,0.0,0.0], 0.21,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
    model.set_value('rip_input', HRU_Input([1.0,0.0,0.0,0.0], 0.02,np.zeros(4), [1,2,3,4], 4, (0,), (0,), 0, np.zeros(4), 0.01, np.zeros(4), 0, 0.0))
    model.set_value('Total_Elevationbands', 4)
    model.set_value('Elevation_Percentage', [0.15,0.26,0.36,0.23])
    model.set_value('bare_storage', Storages(0,np.zeros(4),np.zeros(4),np.zeros(4),0))
    model.set_value('forest_storage', Storages(0,np.zeros(4),np.zeros(4),np.zeros(4),0))
    model.set_value('grass_storage', Storages(0,np.zeros(4),np.zeros(4),np.zeros(4),0))
    model.set_value('rip_storage', Storages(0,np.zeros(4),np.zeros(4),np.zeros(4),0))
  
    # Convert start to date number to be used in if statement
#     start_time = date2num(
#         get_time(CALIBRATION['start']),
#         model.time_units
#     )
    simulated_discharge = []
    # Perform all timesteps of the model, update output fields
    while (model.get_current_time() < (model.get_value_ptr('Date'))[1000]):  
        model.update()
        # Store model time and variable output after the spin up period
        simulated_discharge.append(model.get_value_ptr('Discharge'))
    model.finalize()
    return np.array(simulated_discharge)



def objective_function(parameters, Qobserved, area):
    """Calculate objective function.

    Runs the model, converts the output to GRDC streamflow units
    and calculates NSE from simulation data with observation data.
    This is the function that is going to be optimized by scipy.brute.
    """
    Qmodelled = run_model_cma(parameters)
    QobservedAverage = np.ones(len(Qobserved)) * np.mean(Qobserved)
    Nominator = np.sum((Qmodelled - Qobserved)**2)
    Denominator = np.sum((Qobserved - QobservedAverage)**2)
    NashSutcliffe = 1 - (Nominator / Denominator)
    return 1 - NashSutcliffe


def transform(scaled_parameters):
    """Transforms the scaled_parameter to parameter.

    if x = scaled_parameter and y = parameter,
    then x is in the range [0, 10] and y is in the range [a, b].
    To map the values [0, 10] into [a, b],
    we use the transformations a + (b-a) * x/10.
    For more information on this transformation,
    please see
    http://cma.gforge.inria.fr/cmaes_sourcecode_page.html#practical
    """
    parameters = []
    for scaled_parameter, bound in zip(scaled_parameters, PARAMETERS_BOUNDS):
        scale = (bound[1] - bound[0]) / 10.0
        parameter = bound[0] + (scale * scaled_parameter)
        parameters.append(parameter)
    return parameters


def wrapped_objective_function(scaled_parameters, *args):
    """A wrapper around the objective function.

    The wrapper transforms the scaled_parameters before
    the actual function call.
    """
    parameters = transform(scaled_parameters)
    return objective_function(parameters, *args)


def run_calibration(parameter_bounds, observation, area):
    """Return the optimum parameters found by CMA-ES method."""
    # Set some options for optimization needed when multiprocessing
    options = cma.CMAOptions()
    options.set({
        'bounds': [0, 10],  # for scaled parameters
        'seed': 1234,  # set a seed to reproduce results
        'verbose': -9,  # verbosity of initial/final message: maximally quiet
        'popsize': POPSIZE,
        'maxiter': MAXITER,
        'tolfun': 1e-17,  # very small value due to model creeping behavior
        'tolx': 1e-5
#         'verb_filenameprefix': ('output_dir'),  # cma-es output path
    })

    no_of_variables = len(parameter_bounds)
    # initial mean value and standard deviation
    x0 = 5.0
    sigma0 = 2.0

    # Initialize the CMA-ES
    cma_es = cma.CMAEvolutionStrategy(no_of_variables * [x0], sigma0, options)

    # Store the results of each iteration
    all_scores = []
    all_parameters = []
    
    # Use parallel processing
    with EvalParallel2(number_of_processes=options['popsize']) as evaluations:
        # Use ask-and-tell interface
        while not cma_es.stop():
            solutions = cma_es.ask()
            cma_es.tell(
                solutions,
                evaluations(
                    solutions,
                    fitness_function=wrapped_objective_function,
                    args=(observation, area)
                ),
            )
            # Use transform to return parameters and not scaled ones
            all_parameters.append(transform(cma_es.best.last.x))
            all_scores.append(cma_es.best.last.f)
    print(f"---- CMA-ES stopped due to: {cma_es.stop()} ----")

    # Make full output
    full_output = (transform(cma_es.result.xbest),
                   cma_es.result.fbest,
                   all_parameters,
                   all_scores)
    return full_output


# ### Run calibration

# In[5]:


parameter_bounds = [[-2.0, 2.0],
                     [1.0, 5.0],
                     [0.001, 1.0],
                     [0.1, 0.9],
                     [0.1, 3.0],
                     [0.001, 0.1],
                     [0.4, 0.8],
                     [1.0, 75.0],
                     [0.1, 2.0],
                     [1.0, 10.0],
                     [50.0, 750.0],
                     [0.1, 2.0],
                     [0.1, 5.0],
                     [5.0, 400.0],
                     [0.1, 2.0],
                     [0.1, 8.0],
                     [5.0, 400.0],
                     [0.1, 2.0],
                     [0.1, 3.0],
                     [0.05, 0.5]]

observation = forcing.streamflow[0:1000]
exec_start_time = time.time()
area = 1
full_output = run_calibration(PARAMETERS_BOUNDS, observation, area)
print(f"Calibration is completed in {time.time() - exec_start_time} seconds")





