# -*- coding: utf-8 -*-
"""
Created on Mon May 20 19:02:36 2019

@author: epierre (epierre@triumf.ca)
"""

PSI folder read_me file. v2

The data of the PSI guides experiment can be found on the daq01 computer. It can be accessed using the following command:

ssh -Y ucn@daq01.ucn.triumf.ca

The password is the usual one...

The data are located in the /data/ucn/PSI folder.

This folder contains 3 sub folders:

- data_main contains the data (UCN counts, tof) for the main detector
- data_monitor contains the data (UCN counts, tof) for the monitor detector
- elog_page contains informations about the runs

We will now describe each sub-folder and their dependencies.


data_main folder:

1) the 12 folder is the data tree as it was produced at PSI. The 12 is for the month (data were taken in December). It contains daily folders, named 8, 9, 10, 11 and 12. This 12 folder also contains 3 files: alarmtriggered.txt is an error log, asdf.job is a configuration file for the detector and History.txt contains a summary fo each run taken during the experiment. These files can be ignored for the moment.

Each daily folder contains the data. Each run produces 2 files:
TDDMMYY_XXXX.tof and TDDMMYY_XXXX.txt. DD is for the day the data were taken, MM for the month, YY for the year and XXXX is for the Xth run of the day. For example, the 67th run taken on december 10th, 2017 produced the files T101217_0067.tof and T101217_0067.txt.
The .txt file is a header-summary of the run. It gives indications about the detector (name, serial #, # of pixels), the structure of the .tof file (# of columns, what is each column), # of bins etc. It also summarizes how long the run lasted and how many events were counted.
The .tof files are the data themselves. One file made of X+2 columns (X = # of pixels, it was 16 during this experiment):
the time in bin (1 bin = 0.1 s), the sum of events per time bin (sum of the coming X columns), the # of even per pixel (X columns). If the position doesn't matter, one can just plot the time (bin) as a function of the sum in order to obtain the time of flight spectra.

2) the not_related folder just contains runs which are not related to the experiment. The content of this folder doesn't matter.

3) the sorted_tof_txt contains the dame data as in the 12 folder, but the files are now separated by type (1 folder for the tof files, 1 for the txt files) and not by days. Can be useful for analysis.

4) the sorted_type folder is a work in progress sorting of the data. An efficient way to analyse the data would be to place the relevent files into the right sub_folders. This can be done using the logfiles.

data_monitor folder:

it contains a folder name 12 (like the folder 1)) and all the monotor data (same format as main detector data). Note that the monitor data files were produced at the same time as the main detector data file. Both files have the same name but ar in different folders.

elog_page folder:

it contains a .html and a .pdf version of the logbook. A careful reading of the logbook is needed in order to understand which data are relevant.

