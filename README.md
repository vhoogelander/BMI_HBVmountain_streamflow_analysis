## Climate change impact assessment using a reproducible modeling approach
This repository contains the files and notebooks made for designing a reproducible modeling approach for assessing the impact of climate change on streamflow. The modeling process is designed such that the analysis can be done using three different notebooks: one for forcing generation using the ESMValTool, one for model calibration, and one for running historical and future streamflow simulations and doing the analysis. The modeling process is built on the Julia [HBV-mountain](https://github.com/sarah-hanus/hbv-mountain) hydrological model, for which a Julia [Basic Model Interface](https://github.com/csdms/bmi) (BMI) has been added. This Julia BMI is wrapped with a Python BMI, such that the model can be installed together with [grpc4bmi](https://github.com/eWaterCycle/grpc4bmi) in a Docker container (for later use on the eWaterCycle platform).
Furthermore, a local version of the [ESMValTool](https://github.com/ESMValGroup/ESMValTool) is used, containing ESMValTool recipes specifically for the HBV-mountain model.
The project follows the (FAIR) modeling approaches used in the [eWaterCycle project](https://ewatercycle.readthedocs.io/en/latest/index.html).

## Usage
#### Notebooks
All notebooks can be found in HBVmountain/General. The Example notebook contains an example of how the model can be used for one simulation. The Generate_Forcing notebook is used to generate forcing from the ERA5 or CMIP6 datasets specific for the HBV-mountain model using the ESMValTool version of this repository. The Calibration notebook is used for model calibration and validation, and the Climate_Simulations_Analysis notebook is used for running streamflow simulations and analysis. 

#### Forcing and streamflow observation data
The Data folder contains example data from five different catchments, which is ready for use in the calibration and analysis notebooks. Generating CMIP6 forcings data was done on the HPC cluster Levante from DKRZ, and ERA5 forcings were generated on the eWaterCycle platform (hosted on SURF). For generating forcing data, a catchment shapefile is required (located in the auxiliary_data folder). The model was calibrated using [GRDC](https://www.bafg.de/GRDC/EN/02_srvcs/21_tmsrs/210_prtl/prtl_node.html) streamflow data. 

#### Model set up

To set up the model in a (new) catchment, the following data are needed:
* Catchment shapefile, provided by GRDC together with streamflow data download. 
* [NLCD](https://www.usgs.gov/centers/eros/science/national-land-cover-database) Land cover map (.tiff)
* DEM map (.tiff)
All data must be in WGS84.




Please get in touch for further information.


