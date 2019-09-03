# psi_transmission

This repository contains scripts, Jupyter notebooks, and other documents used for analysis of the UCN guide transmission experiment, carried out in December of 2017 at the [Paul Scherrer Institute](https://www.psi.ch/en). This file provides an overview of its contents.

## data 

Notably missing from the repository are the experimental data, which can be found on the main UCN cluster at the following path: **/ucn/orithyia_data/psi_transmission/name-of-data-folder**

There are three such data folders:

1. **data_ucn** - UCN count data as collected by a main detector and a monitor detector
2. **data_p_beam** - Proton beam current measurements collected throughout the experiment.
2. **data_sim** - Simulation data.

### data_main_detector

There are four important sub-directories to note here: **main_detector/**, **monitor_detector**, **main_detector_sorted/**, and **main_detector_sorted_tof_txt/**.

#### main_detector/

_The following description was provided by Edgard Pierre, a former member of the TUCAN collaboration who led the experimental campaign._

the main_detector folder is the data tree as it was produced at PSI. It contains daily folders, named 8, 9, 10, 11 and 12. This folder also contains 3 files: alarmtriggered.txt is an error log, asdf.job is a configuration file for the detector and History.txt contains a summary fo each run taken during the experiment. These files have been ignored in the analysis thus far.

Each daily folder contains the data. Each run produces 2 files:
TDDMMYY_XXXX.tof and TDDMMYY_XXXX.txt. DD is for the day the data were taken, MM for the month, YY for the year and XXXX is for the Xth run of the day. For example, the 67th run taken on december 10th, 2017 produced the files T101217_0067.tof and T101217_0067.txt.
The .txt file is a header-summary of the run. It gives indications about the detector (name, serial #, # of pixels), the structure of the .tof file (# of columns, what is each column), # of bins etc. It also summarizes how long the run lasted and how many events were counted.
The .tof files are the data themselves. One file made of X+2 columns (X = # of pixels, it was 16 during this experiment):
the time in bin (1 bin = 0.1 s), the sum of events per time bin (sum of the coming X columns), the # of even per pixel (X columns). If the position doesn't matter, one can just plot the time (bin) as a function of the sum in order to obtain the time of flight spectra.

#### monitor_detector/

All of the monitor detector data, in the same format as that of the main detector, but all together in one folder. 

#### main_detector_sorted/

This folder contains both .tof and .txt files for every run that is used in the analysis. This list was produced after inspection of the elog and the data itself. Every file has been renamed according to the convention described in the sorting_data.ipynb notebook file. Use it as a reference.

#### main_detector_sorted_tof_txt/

_also from Edgard..._

the sorted_tof_txt contains the dame data as in the 12 folder, but the files are now separated by type (1 folder for the tof files, 1 for the txt files) and not by days. Can be useful for analysis.

## elog

There are .html, .txt, and .pdf versions of the elog available in the **elog** sub-folder. The naming convention for configurations in the elog is different than my naming convention. The following table presents the conversion:

Elog Name                                     | My Name  |                          
---                                           | ---      | 
first guide (Japanese SS with NiP coating)    | 72 mm SUS Guide with NiP |                           
second guide (Japanese Ti with NiP coating)   | 72 mm SUS Guide with NiP |
SS disk                                       | 72 mm Normalization, No Guide|
two flanges (normally connected to TRIUMF guides) in order to perform a normalization measurement | 85 mm Normalization, no guide |
UGD01                                         | UGD01
UGD03                                         | UGD03
third TRIUMF guide (ep with NiP coating)      | UGD19|
 
This conversion can also be inferred from looking at the run list provided in sorting_data.ipynb

## notebooks and scripts

If trouble is experienced in viewing the notebooks on Github directly, please use [nbviewer](https://nbviewer.jupyter.org/)

**functions.py** file contains every function used in the analysis.

### sorting_data

- a list of runs is provided
- the data file naming conventions are described
- some cuts to individual run data are described

### source_norm

- fits to main detector and monitor detector data are performed with the aim of performing a normalization to the performance of the UCN source

### monitor_detector

- the monitor detector data is analysed

### p_beam_data

- the proton beam current data is analysed

### pre_storage_lifetime

- pre-storage lifetime is analysed

### transmission

- the transmission is calculated for the various guides

### sims

- the simulation data is analysed alongside main detector data

## documents

Some relevant articles, and reports are provided. 

## psi_transmission_report

A git submodule linking to the repository containing a report on the analysis performed. 

## sim_submit

Submission files for simulation. There are config files, materials files, and scripts for running in parallel on the cluster. 
