#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Testing ReadTheDocs

"""
import sys
import os
import numpy as np

from scipy.optimize import curve_fit

import matplotlib as mpl
import matplotlib.pyplot as plt

import time

def func():
    print("Hello world!")
    return 0

# fit function for the degradation of the source performance
def source_fit(t, N_0, A):
    return N_0 + A * t

# fit function for the exponential decay of counts with storage time
def storage_lt_fit(t, N_0, tau):
    return N_0 * np.exp(-t / tau)

# function for getting the start time of the first experimental run of the
# entire campaign
# No Inputs
# Outputs:
#   - the epochal time of the first run returned as an integer
def get_first_run_time():
    
    # grab the txt file for the first run of the experiment,
    #  as it will be used as the t = zero. The monitor runs actually begin 
    # before the main detector runs, so we use the first run of the monitor
    # count data as our time zero.
    f = open('../data_monitor/12/T071217_0001.txt')
    lines = f.readlines()
    f.close()

    # convert start time to epoch time
    date_time = str(7).zfill(2) + '.12.2017 ' + lines[26][15:23]
    pattern = '%d.%m.%Y %H:%M:%S'
    exp_start_time = int(time.mktime(time.strptime(date_time, pattern)))
    return exp_start_time

def load_data(run_list, day_list, normalize_flag=False, just_monitor_flag=False):
    """Loads data from an arbitrary number of runs into numpy arrays
    
    Args:
        run_list (list): list of run numbers (integers)
        day_list (list): list of days (integers) when runs occurred, with 
            one-to-one correspondence to the run_list order.
        normalize_flag (bool, optional): Flag to run normalization routine. 
            Defaults to False.
        just_monitor_flag (bool, optional): Flag to ignore regular data and
            look exclusively at monitor data. Defaults to False.
    
    Returns:
        time_bins (numpy.float64): vector of time data 
        run_times (numpy.float64): vector of run times, these being times in ms
            of the start of each run, in reference to the start of the 
            experiment
        data (numpy.float64): array of main detector count data from all runs
        monitor_data: array of monitor detector count data from all runs
    """

    # !!! Note that I have set the normalize_flag to be false within the 
    # definition until we are satisfied with the normalization routine

    if (normalize_flag):   

        # calculate the UCN yield normalization curve. we can stop doing this 
        # every time data is loaded if speed becomes a problem
        yield_times, yield_result, yield_fit_parameters = UCN_yield()

    # all our data, from every run, has 4000 bins of 1 ms each. We can get away
    # with this sort of lazy hard-coded pre-allocation:
    data = np.zeros((4000, np.size(run_list)))
    monitor_data = np.zeros((4000, np.size(run_list)))

    # first grab the vector of time bins
    time_bins = np.loadtxt("../data_main/12/" + str(day_list[0]) + "/T"
    + str(day_list[0]).zfill(2) + "1217_" + str(run_list[0]).zfill(4) + ".tof", 
    usecols = (0))

    # iterating through the run list
    for i in range(0, np.size(run_list)):
        
        if not (just_monitor_flag):
        
            # grab the vectors of UCN counts
            data[:,i] = np.loadtxt("../data_main/12/" + str(day_list[i]) + "/T"
            + str(day_list[i]).zfill(2) + "1217_" + str(run_list[i]).zfill(4) + 
            ".tof", usecols = (1))

        if (normalize_flag):

            # normalize the vector according to the start time of the run. Here we
            # are normalizing based on the degradation of the UCN source yield over
            # time
            norm_factor = yield_normalization(run_list[i], day_list[i], 
            yield_fit_parameters)

            data[:,i] * norm_factor

        monitor_data[:,i] = np.loadtxt("../data_monitor/12/" + "T"
        + str(day_list[i]).zfill(2) + "1217_" + str(run_list[i]).zfill(4) + 
        ".tof", usecols = (1))

    # get the run times
    run_times = get_run_times(run_list, day_list, just_monitor_flag)   

    return time_bins, run_times, data, monitor_data

config = "NORM"
run_type = "s100"

def load_data_2(config, run_type, normalize_flag = True):
    """A function to load data and sum counts for individual runs.
    
    Arguments:
        config {string} -- A string that determines the experimental 
            configuration of the data to be loaded. The options are:
                'NORM' - normalization
                'JPTI' - JP Ti guide with NiP
                'JPSU' - JP SUS guide with NiP
                'DISK' - SS Disk
                'GD01' - UGD01 guide
                'GD03' - UGD03 guide
                'EPSU' - EP SUS guide with NiP
        run_type {string} -- A string that determines what type of data will be 
            loaded. The options are:
                'shot' - direct shot measurements 
                's005' - 5 second storage
                's020' - 20 second storage
                's100' - 100 second storage
        normalize_flag {boolean, optional} -- Flag to normalize the data 
            accroding to a series of routines dependent on run_type. Defaults
            to True.
    
    Returns:
        numpy.float64 -- An n x 5 data array of the results from loading the
            run data. The number of rows corresponds to the number of runs
            loaded. The three columns are:
                0 - the run start time in seconds since the experimental start
                1 - the storage time (0 if direct shot)
                2 - the number of UCN counts
                3 - sqrt(N) error in number of UCN counts
                4 - [day].[run number] of measurement
    """
    
    # initialize an empty array
    data = np.zeros((1,5))
    
    # perform the source_normalization routine to retrive fit parameters
    if (normalize_flag):
        norm_parameters, norm_errors = source_normalization(run_type)

    # Every file in the directory containing main detector run data is iterated over.
    # Here the /sorted directory contains just runs deemed good for analysis.
    for filename in os.listdir('../data_main/sorted'):

        # Only the files matching our desired configuration and run type are 
        # selected. The '.tof' condition is just so we don't perform the
        # analysis twice per run.
        if ((config in filename) and (run_type in filename) and 
        ('.tof' in filename)):

            f = open( '../data_main/sorted/' + filename[0:22] + '.txt')  
            lines = f.readlines()
            f.close()
            # grab the epoch time for run start
            date_time = filename[1:3].zfill(2) + '.12.2017 ' + lines[26][15:23]
            pattern = '%d.%m.%Y %H:%M:%S'
            run_start_time = int(time.mktime(time.strptime(date_time, pattern)))
            # This function returns the start time of the very first run. 
            run_start_time = run_start_time - get_first_run_time()

            # if the run_type is shot, then storage time is set to 0
            if (run_type == "shot"):

                storage_time = 0

            # otherwise grab the storage time from the run type
            else:

                storage_time = int(run_type[1:4])

            # The data is retrieved from the .tof file
            count_data = np.loadtxt('../data_main/sorted/' + filename[0:22] + 
            '.tof', usecols = (1))

            # specific data cut for run 35 on the 8th
            if ((filename[2:3] == '8') and (filename[10:12] == '35')):

                counts = np.sum(count_data[150:1000])

            # specific data cut for run 66 on the 8th
            elif ((filename[2:3] == '8') and (filename[10:12] == '66')):

                counts = np.sum(count_data[150:1500])

            # specific data cut for run 88 on the 8th
            elif ((filename[2:3] == '8') and (filename[10:12] == '88')):

                counts = np.sum(count_data[150:2500])

            # cut the data normally
            else:

                counts = np.sum(count_data[150:-1])

            # normalize the data depending on the normalize_flag
            if (normalize_flag):
                extrap_counts = source_fit(0, norm_parameters[0], 
                                          norm_parameters[1])
                interp_counts = source_fit(run_start_time, norm_parameters[0], 
                                           norm_parameters[1])
                
                
                norm_factor = extrap_counts / interp_counts
                
                counts = counts * norm_factor
                

            # if this is the first file loaded, then assign the values to the data array, otherwise
            # append the vector of values to the existing array
            if (data[0,0] == 0):

                data[0,0] = run_start_time
                data[0,1] = storage_time
                data[0,2] = counts
                data[0,3] = np.sqrt(counts)

                # saving the [day].[run number] can be useful for debugging
                # requires enlarging the arrays
                data[0,4] = int(filename[1:3]) + \
                (0.001 * int(filename[9:12]))

            # otherwise we make a vector and append it
            else:

                run_data = np.zeros((1,5))
                run_data[0,0] = run_start_time
                run_data[0,1] = storage_time
                run_data[0,2] = counts 
                run_data[0,3] = np.sqrt(counts)

                # saving the [day].[run number] can be useful for debugging
                # requires enlarging the arrays
                run_data[0,4] = int(filename[1:3]) + \
                (0.001 * int(filename[9:12]))
                
                data = np.vstack((data, run_data))
    
    # we return the data sorted by time
    return data[data[:,0].argsort()]

def source_normalization(run_type):
    """Perform a run_type specific fit to data for normalization
    
    Arguments:
        run_type {string} -- A string that determines what type of data will be 
        loaded. The options are:
            'shot' - direct shot measurements 
            's005' - 5 second storage
            's020' - 20 second storage
            's100' - 100 second storage
    
    Returns:
        numpy.float64 -- the computed normalization factor
        numpy.float64 -- error in the normalization factor
    """
        
    # All the runs from the normalization configuration are loaded.
    data = load_data_2('NORM', run_type, normalize_flag = False)

    # The fit is performed.
    popt, pcov = curve_fit(source_fit, data[:,0], data[:,2], p0=[77600, -9], 
    sigma = data[:,3], absolute_sigma = True)

    norm_parameters = popt
    norm_errors     = np.sqrt(np.diag(pcov))

    return norm_parameters, norm_errors

def storage_integrate(data_list):
    """Performs storage time analysis on UCN count data.

    The UCN counts from each given run are integrated, over a set integration
    window dependent on when the valve downstream of the pre-storage volume
    was opened. A fit to determine the pre-storage lifetime is performed.
    
    Args:
        data_list (list): A list of three arrays, each containing run data
        from 100, 20 and 5 second pre-storage time runs. This order matters.
    
    Returns:
        numpy.float64: A 3 x 3 array with rows for each of the three storage
        times, and columns for storage time, UCN counts, and statistical
        uncertainty
    """
    # the result array is 3 x 3.
    result = np.zeros((3,3))

    # the first column is just the storage times
    result[0,0] = 100
    result[1,0] = 20
    result[2,0] = 5

    # !!! the integration range is rather arbitrary, and should be revisited 
    # with a systematic approach? Perhaps based on some fit? For now these
    # hard-coded times will work
    int_times = [110, 35, 20]

    # for each storage time
    for i in range(0, len(data_list)):

        # for each data set (particular run) 
        for j in range(0, np.shape(data_list[i])[1]):

            # sum all the counts
            result[i,1] = result[i,1] + np.sum(data_list[i][int_times[i]:-1,j])

        # divide by the number of runs to get the average
        result[i,1] = result[i,1] / np.shape(data_list[i][1])

        # the statistical uncertainty is calculated
        # !!!i

    # Plotting the result
    plt.scatter(result[:,0], result[:,1], s = 10)
    plt.xlabel('Storage Time [s]')
    plt.ylabel('UCN Counts')
    plt.yscale('log')

    popt, pcov = curve_fit(storage_lt_fit, result[:,0], result[:,1], 
    p0=[150000,28])
    plt.plot(result[:,0], storage_lt_fit(result[:,0], *popt), \
            'r-')
    print("Fit parameters:\n N_0: %.6f,\n tau:   %.6f" % (popt[0], popt[1]))

    # the array returned is sorted by the storage time, for plotting and
    # fitting convenience
    return result

# Function to get the run start times of every run provided. All run times are
# as time elapsed, in seconds, since the first run of the experiment.
# Inputs: 
#   - run numbers as a list
#   - day of runs as an integer (8,9,10,11), in a list with one-to-one
#       correspondence to the run numbers list 
# Outputs:
#   - a vector of times in the order of the provided run number list, the times 
#       are in elapsed seconds since the the start of the first run of the 
#       entire experimental campaign
def get_run_times(run_list, day_list, just_monitor_flag = False):
    
    # initialize an empty array to hold the run start times
    times = np.zeros((len(run_list), 1))
    
    # the start time of the first run is grabbed
    exp_start_time = get_first_run_time()
                     
    for i in range(0, len(run_list)):
        
        # This first if/else case handles the different directory structure for the 
        # monitor v. main detector files
        if (just_monitor_flag):
            
            # read the .txt file for this run
            f = open("../data_monitor/12/" + "T"
            + str(day_list[i]).zfill(2) + "1217_" + str(run_list[i]).zfill(4) + 
            ".txt")         
            
        else:
            # read the .txt file for this run
            f = open("../data_main/12/" + str(day_list[i]) + "/T"
            + str(day_list[i]).zfill(2) + "1217_" + str(run_list[i]).zfill(4) + 
            ".txt")
        
        lines = f.readlines()
        f.close()
        # grab the epoch time for run start
        date_time = str(day_list[i]).zfill(2) + '.12.2017 ' + lines[26][15:23]
        pattern = '%d.%m.%Y %H:%M:%S'
        run_start_time = int(time.mktime(time.strptime(date_time, pattern)))
        
        # record time elapsed since first run in the array
        times[i] = run_start_time - exp_start_time

    return times

# Function to analyze the performance of the UCN source by looking at the direct shot measurements
# performed throughout the experiment.
# Inputs: 
#   - run numbers as a list
#   - day of runs as an integer (8,9,10,11), in a list with one-to-one
#       correspondence to the run numbers list 
# Outputs:
#   - a vector of run times. For each run the run start time is converted to the elapsed number
#     of seconds since the start time of the first run in the run list. 
#   - a vector of integrated counts for each of the runs.
def UCN_yield(plotting_flag = False):
    
    # these are the runs that use the "Two flanges" configuration, which is 
    # returned to intermittently throughout the experimental run. These lists
    # may need to be changed to update this yield calculation routine.
    run_list = [175, 122, 170, 249]
    day_list = [  9,  10,  10,  10]

    # initialize an array to hold the result
    result = np.zeros((len(run_list), 1))
    
    # load data for all the runs into an array
    time_bins, data, monitor_data = load_data(run_list, day_list, False)
    
    for i in range(0, len(run_list)):
        
        # sum all the counts from this direct shot run
        result[i] = np.sum(data[:,i])
        
    # get the run times
    times = get_run_times(run_list, day_list)

    # perform a linear fit
    popt, pcov = curve_fit(source_fit, times[:,0], result[:,0], p0=[77600, -9])

    if (plotting_flag):
        # plotting the result
        plt.scatter(times, result)
        plt.xlabel('Elapsed Time [s]')
        plt.ylabel('Direct Shot UCN Counts')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        plt.plot(times, source_fit(times, *popt), \
        'r-')

    # !!! making this iterative and responsive is a low priority task, can be for
    # fun some time
    # saving the fit parameters, and their respective errors, to an array
    fit_parameters = np.zeros(np.shape(pcov))

    for i in range(0, np.size(popt)):
        fit_parameters[i,0] = popt[i]
        fit_parameters[i,1] = np.sqrt(np.diag(pcov))[i]

    if (plotting_flag):
        # printing the fit parameters and their errors:
        print("Fit parameters:\n N_0 = {} +/- {},\n A: {} +/- {}".format(
                fit_parameters[0,0],
                fit_parameters[0,1],
                fit_parameters[1,0],
                fit_parameters[1,1]))

    # # calculating the error on the normalization factor with standard error
    # # propagation. Recall the source_fit function: return N_0 + A * t
    # norm_factor_err = np.sqrt(N_0_err**2 + A_err**2)

    # print("Error in the normalization factor; to be applied to all runs: %.6f" %
    #     norm_factor_err)

    return times, result, fit_parameters

# Function to calculate a normalization factor based on the degradation of the 
# UCN source yield over time.
# Inputs: 
#   - run_num (int): the number of the run to be normalized
#   - day_num (int): the number of the day of the given run
#   - the fit parameters extracted from fitting the UCN yield over time
# Outputs:
#   - norm_factor (float): the normalization factor to be multiplied by the 
# integrated counts found for the run 
def yield_normalization(run_num, day_num, yield_fit_parameters):
    
    # grab the txt file for the run, as it will be used as the t = zero                 
    f = open("../data_main/12/" + str(day_num) + "/T"
        + str(day_num).zfill(2) + "1217_" + str(run_num).zfill(4) + 
        ".txt")
    lines = f.readlines()
    f.close()
    # convert start time to epoch time
    date_time = str(day_num).zfill(2) + '.12.2017 ' + lines[26][15:23]
    pattern = '%d.%m.%Y %H:%M:%S'
    start_time = int(time.mktime(time.strptime(date_time, pattern)))
    
    # calculate elapsed time since the start of all runs
    # note that this will need to be adjusted according to any updates to the 
    # routine in which UCN yield degradation is analysed.
    exp_start_time = get_first_run_time()
    time_elapsed = start_time - exp_start_time

    # !!!
    norm_factor = source_fit(time_elapsed, yield_fit_parameters[0,0], 
                                yield_fit_parameters[1,0])
    
    return norm_factor

# Function to take 6 lists of runs, 3 for normalization and 3 with a guide and 
# performs a complete transmission analysis
# Inputs: 
#   - run_list_n100, run_list_n20, run_list_n5: lists of normalization runs
#   - run_list_g100, run_list_g20, run_list_g5: lists of guide runs
#   - day_list_n100, day_list_n20, day_list_n5: lists of normalization run days
#   - day_list_g100, day_list_g20, day_list_g5: lists of guide run days
#   - The above 12 lists must be themselves in a list object provided as the 
# single input parameter
#   - The order is crucial. 
# Outputs:
#   - result: 3 x 1 array that contains the 100, 20, and 5 second storage data 
# transmission ratios
def transmission(run_list):
    
    # initialize an array to hold the three transmission results
    result = np.zeros((3,1))
    storage_times = [100, 20, 5]
    
    # !!! the integration range is rather arbitrary, and should be revisited 
    # with a systematic approach? Perhaps based on some fit? For now these
    # hard-coded times will work
    int_times = [110, 35, 20]
    
    # for each list, calculate the average integrated counts after storage
    # compute the ratio between the normalization and guided counts
    # return 3 transmission ratios, one for each of the three storage times. 
    
    for i in range(0, 3):
        
        # variables to hold summed counts
        norm_sum = 0
        guide_sum = 0
    
        # import the normalization data
        time_bins_n, data_n, monitor_data_n = load_data(run_list[i], 
        run_list[i+6])
        
        # import the guide data
        time_bins_g, data_g, monitor_data_g = load_data(run_list[i+3], 
        run_list[i+9])
        
        # we integrate counts in the same fashion as for storage lifetime
        # for each data set (particular run) 
        for j in range(0, np.shape(data_n)[1]):

            # sum all the counts for the normalization and guide configurations, separately
            norm_sum = norm_sum + np.sum(data_n[int_times[i]:-1,j])
            guide_sum = guide_sum + np.sum(data_g[int_times[i]:-1,j])
            
        # divide by the number of runs to get the average
        norm_sum  = norm_sum  / np.shape(data_n)[1]
        guide_sum = guide_sum / np.shape(data_n)[1]
                       
        # log the ratio to the transmission result array
        result[i] = guide_sum / norm_sum

        # the statistical uncertainty is calculated
        # !!!
        
    plt.scatter(storage_times, result * 100)
    plt.xlabel('Storage Time [s]')
    plt.ylabel('Transmission [% / m]')
    
    return storage_times, result

def load_all_data():
    
    norm_data_shot = load_data_2('NORM', 'shot', normalize_flag = True)
    norm_data_5    = load_data_2('NORM', 's005', normalize_flag = True)
    norm_data_20   = load_data_2('NORM', 's020', normalize_flag = True)
    norm_data_100  = load_data_2('NORM', 's100', normalize_flag = True)
    norm_data_list = [norm_data_5, norm_data_20, norm_data_100, norm_data_shot]

    jpsu_data_shot = load_data_2('JPSU', 'shot', normalize_flag = True)
    jpsu_data_5    = load_data_2('JPSU', 's005', normalize_flag = True)
    jpsu_data_20   = load_data_2('JPSU', 's020', normalize_flag = True)
    jpsu_data_100  = load_data_2('JPSU', 's100', normalize_flag = True)
    jpsu_data_list = [jpsu_data_5, jpsu_data_20, jpsu_data_100, jpsu_data_shot]

    jpti_data_shot = load_data_2('JPTI', 'shot', normalize_flag = True)
    jpti_data_5    = load_data_2('JPTI', 's005', normalize_flag = True)
    jpti_data_20   = load_data_2('JPTI', 's020', normalize_flag = True)
    jpti_data_100  = load_data_2('JPTI', 's100', normalize_flag = True)
    jpti_data_list = [jpti_data_5, jpti_data_20, jpti_data_100, jpti_data_shot]

    disk_data_shot = load_data_2('DISK', 'shot', normalize_flag = True)
    disk_data_5    = load_data_2('DISK', 's005', normalize_flag = True)
    disk_data_20   = load_data_2('DISK', 's020', normalize_flag = True)
    disk_data_100  = load_data_2('DISK', 's100', normalize_flag = True)
    disk_data_list = [disk_data_5, disk_data_20, disk_data_100, disk_data_shot]

    gd01_data_shot = load_data_2('GD01', 'shot', normalize_flag = True)
    gd01_data_5    = load_data_2('GD01', 's005', normalize_flag = True)
    gd01_data_20   = load_data_2('GD01', 's020', normalize_flag = True)
    gd01_data_100  = load_data_2('GD01', 's100', normalize_flag = True)
    gd01_data_list = [gd01_data_5, gd01_data_20, gd01_data_100, gd01_data_shot]

    gd03_data_shot = load_data_2('GD03', 'shot', normalize_flag = True)
    gd03_data_5    = load_data_2('GD03', 's005', normalize_flag = True)
    gd03_data_20   = load_data_2('GD03', 's020', normalize_flag = True)
    gd03_data_100  = load_data_2('GD03', 's100', normalize_flag = True)
    gd03_data_list = [gd03_data_5, gd03_data_20, gd03_data_100, gd03_data_shot]

    epsu_data_shot = load_data_2('EPSU', 'shot', normalize_flag = True)
    epsu_data_5    = load_data_2('EPSU', 's005', normalize_flag = True)
    epsu_data_20   = load_data_2('EPSU', 's020', normalize_flag = True)
    epsu_data_100  = load_data_2('EPSU', 's100', normalize_flag = True)
    epsu_data_list = [epsu_data_5, epsu_data_20, epsu_data_100, epsu_data_shot]
    
    return norm_data_list, jpsu_data_list, jpti_data_list, disk_data_list, gd01_data_list, gd03_data_list, epsu_data_list

def storage_lifetime(data_list):
    
    # initialize an array to hold the run averages
    storage_results = np.zeros((3,3))
    
    for i in range(0, 3):
    
        # plot all of the runs together
        plt.errorbar(data_list[i][:,1], data_list[i][:,2], yerr = data_list[i][:,3], fmt = '.', label = '{} s'.format(data_list[i][0,1]))
        
        # compute the standard mean
        storage_results[i,0] = data_list[i][0,1]
        storage_results[i,1] = np.mean(data_list[i][:,2])
        storage_results[i,2] = np.std(data_list[i][:,2])
        
    plt.yscale('log')
    plt.ylabel('UCN Counts');
    plt.xlabel('Storage time [s]')
    plt.legend();
    plt.show()
    
    plt.clf()
    popt, pcov = curve_fit(storage_lt_fit, storage_results[:,0], storage_results[:,1], sigma = storage_results[:,2], p0=[150000, 100], absolute_sigma = True)
    plt.plot(np.linspace(0,100,1000), storage_lt_fit(np.linspace(0,100,1000), *popt));
    plt.errorbar(storage_results[:,0], storage_results[:,1], yerr = storage_results[:,2], fmt = '.')
    plt.ylabel('UCN Counts');
    plt.xlabel('Storage time [s]');
    plt.show()
    plt.clf()
    
    # plotting again with log scale
    plt.plot(np.linspace(0,100,1000), storage_lt_fit(np.linspace(0,100,1000), *popt));
    plt.errorbar(storage_results[:,0], storage_results[:,1], yerr = storage_results[:,2], fmt = '.')
    plt.ylabel('UCN Counts');
    plt.xlabel('Storage time [s]');
    plt.yscale('log')
    
    # printing the fit parameters and their errors:
    fit_parameters = np.zeros(np.shape(pcov))
    for i in range(0, np.size(popt)):
        fit_parameters[i,0] = popt[i]
        fit_parameters[i,1] = np.sqrt(np.diag(pcov))[i]
    print("Fit parameters:\n N_0 = {} +/- {},\n TAU: {} +/- {}".format(
            fit_parameters[0,0],
            fit_parameters[0,1],
            fit_parameters[1,0],
            fit_parameters[1,1]))

    chi_sq_over_dof = np.sum(((storage_results[:,1] - storage_lt_fit(storage_results[:,0], 
                                                                     *popt)) 
                              / storage_results[:,2])**2) / (np.shape(storage_results)[0] - len(popt))
    
    print("chi_sq / dof = {}".format(chi_sq_over_dof))

    
    return storage_results