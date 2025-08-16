# code_CLMU_traffic

## Introduction

This repository is supplementary to the manuscript "Sun, Y., Oleson, K. W., & Zheng, Z. (2025). Modeling urban traffic heat flux in Community Earth System Model".

The objectives of this project are:

- Modify the CESM source code to incorporate an urban traffic heat flux module for quantifying anthropogenic heat flux (AHF);
- Validate model performance with the new traffic module at two sites;
- Examine traffic-induced thermal impacts on the urban environment and building energy use.

## Script and data

### [1_code_modification](./1_code_modification)

The standard source code comes from [CTSM](https://github.com/ESCOMP/CTSM), with the release tag: [ctsm5.3.024](https://github.com/ESCOMP/CTSM/tree/ctsm5.3.024). See modified code lines labeled with `!YS`.

- Create a new module to read time-varying traffic data and provide traffic-related functions:
  - [src/biogeophys/UrbanVehicleType.F90](./1_code_modification/src/biogeophys/UrbanVehicleType.F90)

- Calculate the impervious road width and number of lanes:
  - [src/biogeophys/UrbanParamsType.F90](./1_code_modification/src/biogeophys/UrbanParamsType.F90)

- Pass traffic inputs to compute the traffic heat flux:
  - [src/biogeophys/EnergyFluxType.F90](./1_code_modification/src/biogeophys/EnergyFluxType.F90)
  - [src/biogeophys/UrbanFluxesMod.F90](./1_code_modification/src/biogeophys/UrbanFluxesMod.F90)

- Add the module in the model initialization and computation processes:
  - [src/main/clm_instMod.F90](./1_code_modification/src/main/clm_instMod.F90)
  - [src/main/clm_driver.F90](./1_code_modification/src/main/clm_driver.F90)
- Add time-varying input data to the namelists:
  - [bld/CLMBuildNamelist.pm](./1_code_modification/bld/CLMBuildNamelist.pm)
  - [bld/namelist_files/namelist_defaults_ctsm.xml](./1_code_modification/bld/namelist_files/namelist_defaults_ctsm.xml)
  - [bld/namelist_files/namelist_definition_ctsm.xml](./1_code_modification/bld/namelist_files/namelist_definition_ctsm.xml)

### [2_single_point_simulations](./2_single_point_simulations)

We conducted a pair of single-point simulations (CNTL and TRAF) at the [Capitole of Toulouse, France (FR-Capitole)](./2_single_point_simulations/FR-Capitole) and [Manchester, UK (UK-Manchester)](./2_single_point_simulations/UK-Manchester).

- FR-Capitole
  - [SourceMods](./2_single_point_simulations/FR-Capitole/SourceMods): code used at FR-Capitole simulations with additional modifications to use local parameters. 
  - [datm_files](./2_single_point_simulations/FR-Capitole/datm_files/): atmospheric forcing data, derived from the [Urban-PLUMBER](https://urban-plumber.github.io/).
  - [input_files](./2_single_point_simulations/FR-Capitole/input_files): surface data, derived from the [Urban-PLUMBER](https://urban-plumber.github.io/) and [UTC19 traffic dataset](https://utd19.ethz.ch/).
- UK-Manchester
  - [SourceMods](./2_single_point_simulations/UK-Manchester/SourceMods): code used at UK-Manchester simulations.
  - [datm_files](./2_single_point_simulations/UK-Manchester/datm_files/): atmospheric forcing data from bias-corrected [ERA5-Land reanalysis data](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land).
  - [input_files](./2_single_point_simulations/UK-Manchester/datm_files/): surface data from the [CTSM's default land surface dataset](https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/surfdata_esmf/ctsm5.3.0/) and [Transport for Greater Manchester (TfGM)](https://tfgm.com/)

### [3_simulation_output_analysis](./3_simulation_output_analysis)

The scripts listed below are used to visualize two sites with corresponding traffic diurnal cycles.

| Num. | Subject                                                      | Data process                                                 | Visualization                                                |
| ---- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 3.1  | [Metadata at FR-Capitole](./3_simulation_output_analysis/3.1_FR-Capitole_metadata) | Use [Export.ipynb](./3_simulation_output_analysis/3.1_FR-Capitole_metadata/Export.ipynb) to get the diurnal mean vehicle volume | [Figure.ipynb](./3_simulation_output_analysis/3.1_FR-Capitole_metadata/Figure.ipynb) |
| 3.2  | [Metadata at UK-Manchester](./3_simulation_output_analysis/3.2_UK-Manchester_metadata) | Use [Export.ipynb](/3_simulation_output_analysis/3.2_UK-Manchester_metadata/Export.ipynb) to get the diurnal mean vehicle volume | [Figure.ipynb](./3_simulation_output_analysis/3.2_UK-Manchester_metadata/Figure.ipynb) |

The scripts listed below are used for processing CNTL (`urban_traffic=.false.`) and TRAF (`urban_traffic=.true.`) simulation output and visualization.

| Num. | Subject                                                      | Output data process                                          | Visualization                                                |
| ---- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 3.3  | [Model validation at FR-Capitole](./3_simulation_output_analysis/3.3_FR-Capitole_model_validation/) | Use [Export.ipynb](./3_simulation_output_analysis/3.3_FR-Capitole_model_validation/Export.ipynb) to get monthly mean and hourly mean variables in comparison with observations | [Figure.ipynb](./3_simulation_output_analysis/3.3_FR-Capitole_model_validation/Figure.ipynb) |
| 3.4  | [Anthropoegnic heat at FR-Capitole](./3_simulation_output_analysis/3.4_FR-Capitole_ahf) | Use [Export.ipynb](./3_simulation_output_analysis/3.4_FR-Capitole_ahf/Export.ipynb) to get monthly mean and hourly mean variables related to AHF | [Figure.ipynb](./3_simulation_output_analysis/3.4_FR-Capitole_ahf/Figure.ipynb) |
| 3.5  | [Model validation at UK-Manchester](./3_simulation_output_analysis/3.5_UK-Manchester_model_validation/) | Use [Export.ipynb](./3_simulation_output_analysis/3.5_UK-Manchester_model_validation/Export.ipynb) to get monthly mean and hourly mean variables in comparison with observations | [Figure.ipynb](./3_simulation_output_analysis/3.5_UK-Manchester_model_validation/Figure.ipynb) |
| 3.6  | [Heat stress at UK-Manchester](./3_simulation_output_analysis/3.6_UK-Manchester_heat_stress/) | Use [Export.ipynb](./3_simulation_output_analysis/3.6_UK-Manchester_heat_stress/Export.ipynb) to get heat stress indices | [Figure.ipynb](./3_simulation_output_analysis/3.6_UK-Manchester_heat_stress/Figure.ipynb) |
| 3.7  | [Compare temperatures](./3_simulation_output_analysis/3.7_compare_temperatures/) | Use [Export.ipynb](./3_simulation_output_analysis/3.7_compare_temperatures/Export.ipynb) to get the difference in temperature between TRAF and CNTL simulations | [Figure.ipynb](./3_simulation_output_analysis/3.7_compare_temperatures/Figure.ipynb) |
| 3.8  | [Compare monthly AHF](./3_simulation_output_analysis/3.8_compare_monthly_ahf/) | Use [Export.ipynb](./3_simulation_output_analysis/3.8_compare_monthly_ahf/Export.ipynb) to get the monthly mean AHF from CNTL, TRAF, and [AH4GUC](https://figshare.com/articles/dataset/Global_1-km_present_and_future_hourly_anthropogenic_heat_flux/12612458/6) | [Figure.ipynb](./3_simulation_output_analysis/3.8_compare_monthly_ahf/Figure.ipynb) |
| 3.9  | [Model sensitivity](./3_simulation_output_analysis/3.8_model_sensitivity) | Use [Export.ipynb](./3_simulation_output_analysis/3.8_model_sensitivity/Export.ipynb) to export Taylor Diagram metrics | [Figure.ipynb](./3_simulation_output_analysis/3.8_model_sensitivity/Figure.ipynb) |

### [4_illustration](./4_illustration)

The figures listed below are used to illustrate the details of the model workflow and mechanism.

| Subject                                                      | Visualization                                                |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| Timeline of incorporating AHF in global simulation           | [Figure](./4_illustration/literature_timeline.pdf)           |
| Workflow of incorporating urban traffic modeling in CTSM     | [Figure](./4_illustration/workflow.pdf)                      |
| Biogeophysical mechanism of traffic-induced thermal effects  | [Figure](./4_illustration/mechanism.pdf)                     |
| Traffic-induced changes in heat flux and temperatures at FR-Capitole | [Figure](./4_illustration/FR-Capitole_seasonal_difference.pdf) |
| Traffic-induced changes in heat flux and temperatures at UK-Manchester | [Figure](./4_illustration/UK-Manchester_seasonal_difference.pdf) |
| Community Land Model                                         | [Figure](./4_illustration/CLM.pdf)                           |

### [5_suplimentary_information](./5_suplimentary_information)

The scripts listed below show supplementary information such as input data and simulation results.

| Num. | Subject                                                      | Analysis       | Visualization                                                |
| ---- | ------------------------------------------------------------ | -------------- | ------------------------------------------------------------ |
| 5.1  | [Global number of lanes](./5_suplimentary_information/5.1_global_number_of_lanes/) | Not applicable | [Figure.ipynb](./5_suplimentary_information/5.1_global_number_of_lanes/Figure.ipynb) |

## Acknowledgements

- This work was supported by the Natural Environment Research Council [grant number UKRI1294], and the UKRI Harmonised Impact Acceleration Account, funded via the Economic & Social Research Council (Grant REF: ES/X004759/1) and Engineering & Physical Sciences Research Council (Grant Ref: EP/X525753/1).
- This work used the [ARCHER2 UK National Supercomputing Service](https://www.archer2.ac.uk) and [JASMIN, the UK’s collaborative data analysis environment](https://www.jasmin.ac.uk). 
- We gratefully acknowledge Transport for Greater Manchester (TfGM) for providing traffic data to support this research.
- We appreciate Dr. Xiaodan Xu from Lawrence Berkeley National Laboratory for her valuable insights. 
- [Z. Z.](https://github.com/zhonghua-zheng) appreciates the support provided by the academic start-up funds from the Department of Earth and Environmental Sciences at The University of Manchester. 
- [Y. S.](https://github.com/YuanSun-UoM) is supported by Zhonghua Zheng's academic start-up funds.
- Contributions from [K. W. O.](https://staff.ucar.edu/users/oleson) are based upon work supported by the NSF National Center for Atmospheric Research, which is a major facility sponsored by the U.S. National Science Foundation under Cooperative Agreement No. 1852977.
- The authors declare no conflict of interest.