from BMI_HBVmountain_Python import *
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

def run_model_cma(parameters, forcing, path_to_shapefile, path_to_dem, path_to_nlcd, end_year):

    ob_list = []
    params_list = []
    

   
    model = BMI_HBVmountain(forcing, path_to_shapefile, path_to_dem, path_to_nlcd)
    config_file = model.setup(forcing_netcdf=forcing, bare_parameters=  Parameters(parameters[8], parameters[6], 0, 0, parameters[4],
                                                                     parameters[1], parameters[2], parameters[3], parameters[7], parameters[0]),
                                        forest_parameters=Parameters(parameters[11], parameters[6], 0, parameters[9], parameters[4],
                                                                     parameters[1], parameters[2], parameters[3], parameters[10], parameters[0]),
                                        grass_parameters= Parameters(parameters[14], parameters[6], 0, parameters[12], parameters[4],
                                                                     parameters[1], parameters[2], parameters[3], parameters[13], parameters[0]),
                                        rip_parameters=   Parameters(parameters[17], parameters[6], 0, parameters[15], parameters[18],
                                                                     parameters[1], parameters[2], parameters[3], parameters[16], parameters[0]),
                                        slow_parameters=  Slow_Paramters(parameters[5], parameters[19]))
    model.initialize(config_file)




    Discharge = []
    timestamp = []
    while (model.get_value_ptr('Current_Date') < (datetime.date(end_year, 12, 31))):  
        model.update()
        timestamp.append(model.get_value_ptr('Current_Date'))
        Discharge.append(model.get_value_ptr('Discharge'))

    simulated_discharge_df = pd.DataFrame(
            {'streamflow': Discharge},
            index=pd.to_datetime(timestamp)
        )
    
    precipitation = generate_forcing_from_NETCDF(forcing).prec
    model.finalize()
    return simulated_discharge_df.streamflow, precipitation



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
    PARAMETERS_BOUNDS = [[-2.0, 2.0],
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
    parameters = []
    for scaled_parameter, bound in zip(scaled_parameters, PARAMETERS_BOUNDS):
        scale = (bound[1] - bound[0]) / 10.0
        parameter = bound[0] + (scale * scaled_parameter)
        parameters.append(parameter)
    return parameters

def wrapped_objective_function(scaled_parameters, observation, forcing, path_to_shapefile, path_to_dem, path_to_nlcd,start_year, years_warming_up, end_year):
    """A wrapper around the objective function.

    The wrapper transforms the scaled_parameters before
    the actual function call.
    """
    parameters = transform(scaled_parameters)
    Qmodelled, precipitation = run_model_cma(parameters, forcing, path_to_shapefile, path_to_dem, path_to_nlcd, end_year)
    precipitation.index, Qmodelled.index, observation.index = pd.to_datetime(precipitation.index), pd.to_datetime(Qmodelled.index), pd.to_datetime(observation.index)
    mask = (observation.index >= Qmodelled.index[0]) & (observation.index <= Qmodelled.index[-1])
    observation = observation.loc[mask]

    ED = multi_objective(Qmodelled.loc[Qmodelled.index.year >= (start_year+years_warming_up)], observation.loc[observation.index.year >= (start_year+years_warming_up)], precipitation)[0] 
    return ED

def cma_calibration(parameter_bounds, observation, MAXITER, forcing, path_to_shapefile, path_to_dem, path_to_nlcd, start_year, years_warming_up, end_year):
    """Return the optimum parameters found by CMA-ES method."""
    # Set some options for optimization needed when multiprocessing
    options = cma.CMAOptions()
    POPSIZE = 9
    options.set({
        'bounds': [0, 10],  # for scaled parameters
#         'seed': 1234,  # set a seed to reproduce results
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
    # with EvalParallel2(number_of_processes=options['popsize']) as evaluations:
        # Use ask-and-tell interface
    while not cma_es.stop():
        solutions = cma_es.ask()
        cma_es.tell(
            solutions,

            [wrapped_objective_function(x, observation, forcing, path_to_shapefile, path_to_dem, path_to_nlcd, start_year, years_warming_up, end_year) for x in solutions]
                #function_values=wrapped_objective_function(solutions, observation, area)
            )#,

        # Use transform to return parameters and not scaled ones
        all_parameters.append(transform(cma_es.best.last.x))
        all_scores.append(cma_es.best.last.f)
#     print(f"---- CMA-ES stopped due to: {cma_es.stop()} ----")

    # Make full output
    full_output = (transform(cma_es.result.xbest),
                   cma_es.result.fbest,
                   all_parameters,
                   all_scores)
    return full_output


def full_calibration(ntimes, observation, MAXITER, forcing, path_to_shapefile, path_to_dem, path_to_nlcd, start_year, years_warming_up, end_year):
    exec_start_time = time.time()
    PARAMETERS_BOUNDS = [[-2.0, 2.0],
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
    area = gpd.read_file(path_to_shapefile).area_hys.values[0] #km2
    observation = pd.read_csv(observation, index_col=0).streamflow / (area * 1e6) * 1000 *86400
    ob_list = []
    params_list = []
    results_doc_list = []
    for n in range(ntimes):
        full_output = cma_calibration(PARAMETERS_BOUNDS, observation, MAXITER, forcing, path_to_shapefile, path_to_dem, path_to_nlcd, start_year, years_warming_up, end_year)
        best_par, best_score, pars, scores = full_output
        
        simulated_discharge, precipitation = run_model_cma(best_par, forcing, path_to_shapefile, path_to_dem, path_to_nlcd, end_year)
        mask = (observation.index >= simulated_discharge.index[0]) & (observation.index <= simulated_discharge.index[-1])
        observation = observation.loc[mask]
        
        objective_function = multi_objective(simulated_discharge.loc[simulated_discharge.index.year >= (start_year+ years_warming_up)], observation.loc[observation.index.year >= (start_year+ years_warming_up)], precipitation)
        ob_list.append(objective_function)
        params_list.append(best_par)   
        

        
        
        other_values = [
            ([float(z) for z in x], float(1 - y)) for x, y in zip(pars, scores)
        ]
        results_doc = {
            'title': f'HBVmountain model optimized parameter based on ED,{n+1}st simulation',
            'best_parameter_value': [float(x) for x in best_par],
            'best_NSE_value': float(1 - best_score),
            'other_values': other_values,
        }

        results_doc_list.append(results_doc)
    paramset = pd.DataFrame(columns=['ED', 'NSE', 'logNSE', 'NSEfdc', 'NSErunoff', 'Temp_Thresh', 'Meltfactor', 'Mm', 'Ratio_Pref', 'Kf', 'Ks', 'Ce', 'Soilstoragecapacity_Bare', 'beta_Bare', 'Interceptioncapacity_Forest', 'Soilstoragecapacity_Forest', 'beta_Forest', 'Interceptioncapacity_Grass', 'Soilstoragecapacity_Grass', 'beta_Grass', 'Interceptioncapacity_Rip', 'Soilstoragecapacity_Rip', 'beta_Rip', 'Kf_Rip', 'Ratio_Riparian'])
    for i in range(len(ob_list)):
        paramset.loc[i] = [ob_list[i][0], ob_list[i][1], ob_list[i][2], ob_list[i][3], ob_list[i][4], params_list[i][0], params_list[i][1], params_list[i][2], params_list[i][3], params_list[i][4], params_list[i][5], params_list[i][6], params_list[i][7], params_list[i][8], params_list[i][9], params_list[i][10], params_list[i][11], params_list[i][12], params_list[i][13], params_list[i][14], params_list[i][15], params_list[i][16], params_list[i][17], params_list[i][18], params_list[i][19]]
    paramset.ED = 1 - paramset.ED
    filename = 'cma_calibration_thundercreek.yml'
    yaml = YAML()
    yaml.default_flow_style = False
    with open(filename, 'w+') as f:
        yaml.dump(results_doc_list, f)
    print(f'Calibration results saved in {filename}')
    print(f"Calibration is completed in {time.time() - exec_start_time} seconds")
    return paramset




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
    paramset_item = full_calibration(partition_item,  observation, MAXITER, forcing, path_to_shapefile, path_to_dem, path_to_nlcd, start_year, years_warming_up, end_year)
    paramset = comm.gather(paramset_item, root=0)

    if rank == 0:
        paramset = pd.concat(paramset)
        name = 'paramsets_thundercreek.csv'
        outdir = './output'
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        fullname = os.path.join(outdir, name)
        paramset.to_csv(fullname)
   


if __name__ == '__main__':
    exec_start_time = time.time()
    main()
    print(f"Calibration is completed in {time.time() - exec_start_time} seconds")

