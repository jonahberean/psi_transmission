# -*- coding: utf-8 -*-
"""Functions for analysis of 2017 PSI transmission experiment

This module contains a variety of functions that are useful for loading
and analyzing the data set from the 2017 PSI UCN transmission experiment.
"""
# Some standard import statements

import sys
import os
import logging
import time
import decimal
import uncertainties
from uncertainties import unumpy
from uncertainties import *
import numpy as np
from scipy.optimize import curve_fit
from IPython import get_ipython
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

###############################################################################
###############################################################################

def linear_fit(t, N_0, y):
    """function for a linear fit
    
    Arguments:
        t {float} -- time of run
        N_0 {float} -- initial counts at t = 0
        y {float} -- slope of the loss
    
    Returns:
        float -- counts as a function of time; N(t)
    """    
    return N_0 - y * t

###############################################################################
###############################################################################
def storage_lt_fit(t, N_0, tau):
    """function for a storage lifetime fit
    
    Arguments:
        t {float} -- storage time
        N_0 {float} -- initial counts at t = 0
        tau {float} -- characteristic pre-storage decay time
    
    Returns:
        float -- counts as a function of time; N(t)
    """   
    return N_0 * np.exp(-t / tau)

###############################################################################
###############################################################################


def get_start_time():
    """Retrieves the UNIX epoch time stamp of the first proton beam current
    measurement, of the 2 second period data. This is used as a reference, 
    t = 0, time for all analysis
    
    Returns:
        int: the UNIX epoch time stamp
    """
    
    # read the 8th of December data as a list of strings
    f = open('../data_p_beam/2_second/20171208.csv')
    lines = f.readlines()
    f.close()
    
    # convert the measurement time to epoch time
    date_time = lines[1][0:10] + ' ' + lines[1][11:19]
    # print(date_time)
    pattern = '%Y-%m-%d %H:%M:%S'
    start_time = int(time.mktime(time.strptime(date_time, pattern)))
    
    return start_time

###############################################################################
###############################################################################

def load_main(config, run_type):
    """A function to load data and sum counts for runs of a given
        configuration and pre-storage time
    
    Arguments:
        config {string} -- A string that determines the experimental 
            configuration of the data to be loaded. The options are:
                'NOMI' - nominal, guide-less
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
    
    Returns:
        numpy.float64 -- An n x 5 data array of the results from loading the
            run data. The number of rows corresponds to the number of runs
            loaded. The five columns are:
                0 - the run start time in seconds since the experimental start
                1 - the storage time (0 if direct shot)
                2 - the number of UCN counts
                3 - sqrt(N) error in number of UCN counts
                4 - [day].[run number] of measurement
    """
    # start_time is hard-coded here as the UNIX time stamp of the first 
    # proton beam current measurement, of the 2 second data. 
    start_time = get_start_time()

    # instantiate a new numpy array 
    all_data = np.empty((0,5), float)
            
    for filename in os.listdir('../data_ucn/main_detector_sorted'):

        # Only the files matching our desired configuration and run 
        # type are selected. The '.tof' condition is just so we 
        # don't perform the analysis twice per run (since it would
        # otherwise match to the .tof and the .txt files)
        if ((config in filename) and (run_type in filename) and 
        ('.tof' in filename)):

            # grab from the text file associated with the run
            f = open( '../data_ucn/main_detector_sorted/' + filename[0:22] 
                     + '.txt')  
            lines = f.readlines()
            f.close()

            # grab the epoch time for run start
            date_time = filename[1:3].zfill(2) + \
                        '.12.2017 ' + \
                        lines[26][15:23]

            pattern = '%d.%m.%Y %H:%M:%S'
            run_time = int(time.mktime(time.strptime(date_time, pattern)))

            # reset the run_start_time with reference to the
            # t = 0 time
            run_time = run_time - start_time

            # grab the storage time
            if (run_type == 'shot'):

                storage_time = 0

            else:    

                storage_time = int(run_type[1:4])

            # The data is retrieved from the .tof file
            count_data = np.loadtxt('../data_ucn/main_detector_sorted/' + 
                                    filename[0:22] + '.tof',
                                    usecols = (1))

            # various cuts to the time-of-flight spectra are made
            N = spectrum_cuts(filename, count_data, run_type)

            # saving the [day].[run number] can be useful for debugging
            day_run_no = int(filename[1:3]) + (0.001
                                               * int(filename[9:12]))

            # append the loaded data to the existing array
            all_data = np.append(all_data, 
                                [[run_time,
                                   storage_time,
                                   N,
                                   np.sqrt(N),
                                   day_run_no]],
                                axis=0)

    # return the completed array, sorted along the time axis
    return all_data[all_data[:,0].argsort()]

###############################################################################
###############################################################################

def spectrum_cuts(filename, count_data, run_type):
    """makes specific cuts to the time-of-flight spectrum based on known issues
    in the data
    
    Arguments:
        filename {string} -- name of the file being loaded
        count_data {numpy.float64} -- summed counts
    
    Returns:
        numpy.float64 -- summed counts after the cut
    """

    # this if/else sequence handles cuts of the data, which for 
    # some runs is specific based on the experimental 
    # conditions
    # specific data cut for run 35 on the 8th
    if ((filename[2:3] == '8') and 
        (filename[10:12] == '35')):

        N = np.sum(count_data[150:1000])

    # specific data cut for run 66 on the 8th
    elif ((filename[2:3] == '8') and 
            (filename[10:12] == '66')):

        N = np.sum(count_data[150:1500])

    # specific data cut for run 88 on the 8th
    elif ((filename[2:3] == '8') and 
            (filename[10:12] == '88')):

        N = np.sum(count_data[150:2500])

    # if it is a shot run we take all the counts
    elif (run_type == 'shot'):
    
        N = np.sum(count_data)
    
    # otherwise cut the data normally for a pre-storage run
    # this cuts out the initial background appearing from
    # irradiation
    else:

        N = np.sum(count_data[150:-1])

    return N

###############################################################################
###############################################################################

def load_all_main(norm_flag = True):
    """A function to load data and sum counts for all the run data available
    
    Arguments:
        norm_flag {boolean, optional} -- flag to turn on normalization based
        on sD2 losses. Defaults to True.
    
    Returns:
        dict -- a dictionary of arrays of the same structure as is returned by
        load_data(). The key pairs to be used are:
        key 0: config {string} -- The options are:
            'NOMI' - nominal, guide-less
            'JPTI' - JP Ti guide with NiP
            'JPSU' - JP SUS guide with NiP
            'DISK' - SS Disk
            'GD01' - UGD01 guide
            'GD03' - UGD03 guide
            'EPSU' - EP SUS guide with NiP
            'all'  - all of the above
        key 1: run_type {string} -- The options are:
            'shot' - direct shot measurements 
            's005' - 5 second storage
            's020' - 20 second storage
            's100' - 100 second storage
            'all'  - all of the above

        dict -- dictionary of values (unumpy ufloat object) 
        of the results from the sD2 losses normalization. The key pairs to be 
        used are:
        key 0: run_type {string} -- The options are:
            'shot' - direct shot measurements 
            's005' - 5 second storage
            's020' - 20 second storage
            's100' - 100 second storage
        key 1: parameter {string} -- The options are:
            'N_0'     - counts at time 0 +/- error 
            'y'       - loss rate +/- error
    """
    # instantiate configuration and run type lists
    config_list = ['JPTI', 'JPSU', 'DISK', 'GD01', 'GD03', 'EPSU']
    run_type_list = ['shot', 's005', 's020', 's100']

    # instantiate dictionary to hold all main detector data
    data_dict = {}

    # load data for normalization
    if (norm_flag == True):

        # initialize dictionary to hold parameter results
        norm_dict = {}
        
        # iterate over each run type
        for run_type in run_type_list:

            # load the main detector data for the TRIUMF-style normalization
            # configuration
            arr = load_main('NOMI', run_type)

            # get the normaliization fit parameters
            popt, pcov = curve_fit(linear_fit, arr[:,0], arr[:,2], 
                                sigma = arr[:,3], absolute_sigma = True)
            
            # saving the fit results to the dictionary, as uncertainty objects
            norm_dict[run_type, 'N_0'] = ufloat(popt[0], 
                                                np.sqrt(np.diag(pcov))[0])
            norm_dict[run_type, 'y']   = ufloat(popt[1], 
                                                np.sqrt(np.diag(pcov))[1])

            # normalize the very data used to calculate the normalization and
            arr = sD2_normalize(arr, norm_dict, run_type)

            data_dict['NOMI', run_type] = arr

    else:

        norm_dict = None

    # the 'all', run_type dicts and the 'all', 'all' dict must be instantiated
    # prior to the main for loops
    for run_type in run_type_list:

        data_dict['all', run_type] = np.empty((0,5), float)

    data_dict['all', 'all'] = np.empty((0,5), float)

    for config in config_list:

        # at start of each config, initialize the empty 'config', 'all' array
        data_dict[config, 'all'] = np.empty((0,5), float)

        for run_type in run_type_list:

            # load the appropriate data into an array
            arr = load_main(config, run_type)

            # perform the normalization for sD2 losses
            if (norm_flag == True): 

                arr = sD2_normalize(arr, norm_dict, run_type)

            data_dict[config, run_type] = arr

            data_dict['all', run_type] = np.append(data_dict['all', run_type],
                                                arr,
                                                axis = 0)

            data_dict[config, 'all'] = np.append(data_dict[config, 'all'],
                                                arr,
                                                axis = 0)

            data_dict['all', 'all'] = np.append(data_dict['all', 'all'],
                                                arr,
                                                axis = 0)

    return data_dict, norm_dict

###############################################################################
###############################################################################

def load_monitor():
    """A function to load and sum all UCN count data from monitor detector
    runs
    
    Returns:
        numpy.float64 -- an n x 4 array of the resulting data The five columns 
        are:
            0 - the run start time in seconds since the experimental start
            1 - the number of UCN counts
            2 - the Poisson error in UCN counts, \sqrt{N}
            3 - [day].[run number] of measurement
    """    

    # get the start time
    start_time = get_start_time()

    # initialize an array to hold the data
    monitor_data = np.empty((0,4), float)

     # loop through the files and load the data
    for filename in os.listdir('../data_ucn/monitor_detector'):
        
        # get the time stamp from the txt file and the counts from the tof file
        # but we only check for one, so that we don't do each twice.
        if(filename[0] == 'T' and 'tof' in filename):
            
            # print(filename[0:12])

            # grab from the text file associated with the run
            f = open('../data_ucn/monitor_detector/' 
                            + filename[0:12] 
                            + '.txt')  

            lines = f.readlines()
            f.close()

            # grab the epoch time for run start
            date_time = filename[1:3].zfill(2) + '.12.2017 '\
                + lines[26][15:23]
            
            pattern = '%d.%m.%Y %H:%M:%S'
            run_time = int(time.mktime(
                time.strptime(date_time, pattern)))

            # reset the run_start_time with reference to the
            # t = 0 time
            # !!! temporarily use the raw UNIX epoch time stamp
            # run_time = run_time - start_time

            # load the monitor count data
            arr = np.loadtxt('../data_ucn/monitor_detector/' + filename,
                                usecols = (1))

            # sum the counts
            counts = np.sum(arr)

            # saving the [day].[run number] can be useful for debugging
            day_run_no = int(filename[1:3]) + (0.001
                                               * int(filename[9:12]))

            # the current data is appended to the existing data array
            monitor_data = np.append(monitor_data, [[run_time, 
                                                    counts, 
                                                    np.sqrt(counts),
                                                    day_run_no]], axis = 0)
    
    return monitor_data

###############################################################################
###############################################################################

def load_p_beam_10s():
    """Loads all of the 10 second period proton beam data into an array
    
    Returns:
        p_beam_data (numpy.float64): a (2 x n) array of all of the proton beam 
            monitoring data. 
                - row 0: time elapsed in seconds since the first measurement
                - row 1: beam current in uA - monitoring data
                - row 3: beam current in uA - timing data (not to be trusted
                    for absolute value)
    """
    
    # get start time
    start_time = get_start_time()

    # instantiate array to hold the resulting data, empty and single column 
    # at first, for data to be successively stacked
    p_beam_data = np.empty((0,3), float)
    
    # loop through the files and load the data
    for filename in os.listdir('../data_p_beam/10_second'):
        
        # all of the csv file is converted to a list of strings for extracting
        # the time data
        f = open('../data_p_beam/10_second/' + filename)
        lines = f.readlines()
        f.close()
        
        # instantiate an array to hold the measurement times
        arr = np.zeros((np.shape(lines)[0] - 1, 3))
        
        # loop over every row in the csv file, skipping line 1
        for i in range(0, np.shape(arr)[0]):
            
            # convert the measurement time to epoch time
            date_time = lines[i + 1][0:10] + ' ' + lines[i + 1][11:19]
            # print(date_time)
            pattern = '%d.%m.%Y %H:%M:%S'
            measurement_time = int(
                time.mktime(time.strptime(date_time, pattern)))
            
            # save the elapsed time to the times array
            arr[i, 0] = measurement_time - start_time

        # the current data is loaded into a numpy array
        arr[:,1:3] = np.loadtxt('../data_p_beam/10_second/' + filename, 
                              delimiter = ';', 
                              skiprows=1, 
                              usecols=(86,88));
        
        # removing the 0 values
        for i in range(0,np.shape(arr)[0]):

            if (arr[i,1] == 0):

                arr[i,1] = float('nan')

            if (arr[i,2] == 0):

                arr[i,2] = float('nan')


        # append the time and count data to the array
        p_beam_data = np.append(p_beam_data, 
                                arr, axis = 0)

    return p_beam_data

###############################################################################
###############################################################################

def load_p_beam_2s():
    """Loads all of the 2 second period proton beam data into a numpy array
    
    Returns:
        p_beam_data (numpy.float64): a (2 x n) array of all of the proton beam 
            monitoring data. 
                - row 0: time elapsed in seconds since the first measurement
                - row 1: beam current in uA - monitoring data
    """
    # get start time
    start_time = get_start_time()
    
    # instantiate array to hold the resulting data, empty and single column 
    # at first, for data to be successively stacked
    p_beam_data = np.empty((0,2), float)
    
    # loop through the files and load the data
    for filename in os.listdir('../data_p_beam/2_second'):
        
        # all of the csv file is converted to a list of strings for extracting
        # the time data
        f = open('../data_p_beam/2_second/' + filename)
        lines = f.readlines()
        f.close()
        
        # instantiate an array to hold the measurement times
        arr = np.zeros((np.shape(lines)[0] - 1, 2))
        
        # loop over every row in the csv file, skipping line 1
        for i in range(0, np.shape(arr)[0]):
            
            # convert the measurement time to epoch time
            date_time = lines[i + 1][0:10] + ' ' + lines[i + 1][11:19]
            # print(date_time)
            pattern = '%Y-%m-%d %H:%M:%S'
            measurement_time = int(
                time.mktime(time.strptime(date_time, pattern)))
            
            # save the elapsed time to the times array
            arr[i, 0] = measurement_time - start_time

        # the current data is loaded into a numpy array
        arr[:,1] = np.loadtxt('../data_p_beam/2_second/' + filename, 
                              delimiter = ',', 
                              skiprows=1, 
                              usecols=(1));
        
        # removing the 0 values
        for i in range(0,np.shape(arr)[0]):

            if (arr[i,1] == 0):

                arr[i,1] = float('nan')


        # append the time and count data to the array
        p_beam_data = np.append(p_beam_data, 
                                arr, axis = 0)

    return p_beam_data

###############################################################################
###############################################################################

def find_coincidences(p_beam_data, main_data_dict, window, plotting_flag = False):
    """Searches for coincident runs and measurements, from the main detector
        data and the proton beam current data, respectively.
    
    Arguments:
        p_beam_data (numpy.float64) -- a (2 x n) array of all of the proton beam 
            monitoring data. 
                - row 0: time elapsed in seconds since the first measurement
                - row 1: beam current in uA - monitoring data
        main_data_dict (dict) -- see docstring of load_all_data()
        window (int) -- number of seconds for coincidence time window. MUST BE 
            < 400
        plotting_flag (boolean, optional) -- flag to enable plotting. Defaults
            to False.
    
    Returns:
        (dict) -- a dictionary of arrays for coincident proton beam
            measurements, and main neutron count data, 
            labelled by run type. The key options are:
            key 0: data_type {string} -- The options are:
                'proton' - proton beam current measurements
                'neutron' - neutron count data          
            key 1: run_type {string} -- The options are:
                'shot' - direct shot measurements 
                's005' - 5 second storage
                's020' - 20 second storage
                's100' - 100 second storage
            
    """

    # instantiate a new dictionary
    coincidence_dict = {}

    # iterate over run types
    run_type_list = ['shot', 's005', 's020', 's100']

    for run_type in run_type_list:

        # grab the neutron data
        n_arr = main_data_dict['all', run_type]
        
        # grab and sort the proton beam data
        p_arr = p_beam_data
        p_arr = p_arr[p_arr[:,0].argsort()]
        
        # arrays to keep track of indices
        p_indices = np.empty((1,0), float)
        n_indices = np.empty((1,0), float)

        # arrays to hold the coincident points found
        p_coincident = np.empty((0,2), float)
        n_coincident = np.empty((0,5), float)

        # iterating over every time stamp of the proton beam current 
        # measurements
        for i in range(0, np.size(p_arr[:,0])):

            # only where the proton beam current has non-nan values
            if np.isnan(p_arr[i,1]) == False:
                
                # get the index, in the neutron count data, that matches this
                # beam current measurement time stamp
                n_indices_to_add = np.argwhere(np.abs(n_arr[:,0] 
                                                - p_arr[i,0]) < window)
                
                # if its not empty, then a coincidence was found
                if (np.size(n_indices_to_add) != 0):

                    # keep track of indices
                    p_indices = np.append(p_indices, i)
                    n_indices = np.append(n_indices, n_indices_to_add)
                    
                    # add data to coincidence arrays
                    p_coincident = np.append(p_coincident, 
                                            [p_arr[i,:]], axis = 0)
                    n_coincident = np.append(n_coincident, 
                                            n_arr[n_indices_to_add[0], :], 
                                            axis = 0)
                                    
        # save the completed arrays to the coincidence dict
        coincidence_dict['proton', run_type] = p_coincident
        coincidence_dict['neutron', run_type] = n_coincident
                    
        if (plotting_flag):

                # instantiate the subplots 
                fig, ax1 = plt.subplots()

                if(run_type == 'shot'):

                    ax1.set_title('direct shot')

                else:

                    ax1.set_title(run_type[1:4] + ' second storage')

                # plot the proton beam current data
                ax1.scatter(p_coincident[:,0], p_coincident[:,1], color = 'r')

                # presentation stuff
                ax1.set_xlabel('Time Elapsed (s)')
                ax1.set_ylabel(r'Proton Beam Current [$\mu$A]', 
                    color = 'r')
                ax1.tick_params(axis='y')
                ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

                # instantiate a second axes that shares the same x-axis
                ax2 = ax1.twinx()  
                
                # plot the neutron count data for runs of this pre-storage
                # time
                ax2.errorbar(n_arr[:,0], n_arr[:,2], 
                                yerr = n_arr[:,3], 
                                marker = 'o',
                                fmt = '.',
                                color = 'b', 
                                label = 'not coincident')

                # plotting the coincident neutrond data afterwards in black
                # to highlight these runs vs non-coincident runs
                ax2.errorbar(n_coincident[:,0], n_coincident[:,2], 
                                yerr = n_coincident[:,3],
                                marker = 'o',
                                fmt = '.', 
                                color = 'k',
                                label = 'coincident')

                # presentation stuff
                ax2.set_ylabel('UCN Counts')  
                ax2.set_yscale('log')
                ax2.legend()
                fig.tight_layout()  

    return coincidence_dict

###############################################################################
###############################################################################

def sD2_normalize(arr, norm_dict, run_type):
    """normalizes the given array with the norm_dict values provided
    
    Arguments:
        arr {numpy.float64} -- see load_main_all for structure
        norm_dict {dict} -- see load_main_all for key-value details
        run_type {string} -- see load_main_all for options
    
    Returns:
        numpy.float64 -- the normalized array with same structure
    """

    # generate uncertainty array of the run data
    uarr = unumpy.umatrix(arr[:,2], arr[:,3])

    # generate uncertainty float object run time, which has no uncertainty
    run_time = unumpy.umatrix(arr[:,0], 0)

    # interpolate along the normalization fit for each run start time
    uinterp = linear_fit(run_time, 
                    norm_dict[run_type, 'N_0'], 
                    norm_dict[run_type, 'y'])

    # get the normalization factor required for each point
    unorm = norm_dict[run_type, 'N_0'] / uinterp

    # normalize the run data
    uarr = np.multiply(uarr, unorm)

    # update the array values, and hence the dictionary by mutablity in python
    arr[:,2] = unumpy.nominal_values(uarr)
    arr[:,3] = unumpy.std_devs(uarr)

    return arr