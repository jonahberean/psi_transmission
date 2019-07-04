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