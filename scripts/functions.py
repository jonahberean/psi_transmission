"""
!!! complete the docstring
Our first Python module. This initial string is a module-level documentation string.
It is not a necessary component of the module. It is a useful way to describe the
purpose of your module.
"""

# Some standard import statements

import sys
import os
import logging
import time
import decimal
import numpy as np
from scipy.optimize import curve_fit
from IPython import get_ipython
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

def hello():
    print('hello')

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

def load_main(config, run_type, norm_dict_in = None):
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
        norm_dict_in {dict} -- dictionary of values of the results from the
            ucn yield analysis. Defaults to None which avoids normalization.
            The key pairs to be used are:
                key 0: run_type {string} -- The options are:
                'shot' - direct shot measurements 
                's005' - 5 second storage
                's020' - 20 second storage
                's100' - 100 second storage
                key 1: parameter {string} -- The options are:
                'N_0'     - counts at time 0
                'y'       - loss rate
                'N_0_err' - associated error
                'y_err'   - associated error
    
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

    if norm_dict_in != None:
        
        norm_dict = dict(norm_dict_in)

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
            date_time = filename[1:3].zfill(2) + '.12.2017 '\
                + lines[26][15:23]
            pattern = '%d.%m.%Y %H:%M:%S'
            run_time = int(time.mktime(
                time.strptime(date_time, pattern)))

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

            # this if/else sequence handles cuts of the data, which for 
            # some runs is specific based on the experimental 
            # conditions
            # !!! These need to be re-documented somewhere in ipynb
            # !!! The eLog parsing needs to be summarized as well
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
            # !!! But does this mean that we're counting that initial proton 
            # irradiation background in here? What is that background from?
            # should it be chracterized and removed? Does it fluctuate over
            # time? If so could it be tied to the current measurements and then
            # corrected out?
            elif (run_type == 'shot'):
            
                N = np.sum(count_data)
            
            # otherwise cut the data normally for a pre-storage run
            # this cuts out the initial background appearing from
            # irradiation
            else:

                N = np.sum(count_data[150:-1])

            # normalize the data depending on the normalize_flag
            if (norm_dict_in != None):

                # !!! the error associated with the normalization 
                # routine has not been implemented

                # the fits to the nominal configuration data provide
                # a benchmark for percentage loss of absolute counts
                # depending on the time having elapsed since the start
                # of the experiment. 

                denom = linear_fit(run_time, norm_dict[run_type, 'N_0'], 
                                           norm_dict[run_type, 'y'])

                S = norm_dict[run_type, 'N_0'] / denom

                # normalize the counts
                N = N * S

                # compute the uncertainty in the denominator for S, i.e. the 
                # N(t) calculation from the nominal configuration data fit
                N_0_err = norm_dict[run_type, 'N_0_err']
                y_err = norm_dict[run_type, 'y_err']
                denom_err = np.sqrt((N_0_err)**2 + (y_err)**2)

                # calculate the fractional uncertainty in the numerator for S
                N_0_frac_err = N_0_err / norm_dict[run_type, 'N_0']

                # calculate the fractional uncertainty in S
                S_frac_err = np.sqrt((N_0_frac_err)**2 +
                                     (denom_err / denom)**2)

                # calculate the resulting fractional uncertainty in counts

                N_frac_err = np.sqrt((np.sqrt(N) / N)**2 + 
                                                            (S_frac_err)**2)

                # the absolute unceratainty 
                N_err = N_frac_err * N

            # if no normalization, then N_err is just sqrt(N) of Poisson
            else:

                N_err = np.sqrt(N)

            # saving the [day].[run number] can be useful for debugging
            day_run_no = int(filename[1:3]) + (0.001
                                               * int(filename[9:12]))

            # append the loaded data to the existing array
            all_data = np.append(all_data, 
                                [[run_time,
                                   storage_time,
                                   N,
                                   N_err,
                                   day_run_no]],
                                axis=0)

    # return the completed array, sorted along the time axis
    return all_data[all_data[:,0].argsort()]

###############################################################################
###############################################################################

def load_all_main(norm_dict = None):
    """A function to load data and sum counts for all the run data available
    
    Arguments:
        norm_dict {dict} -- dictionary of values of the results from the
            ucn yield analysis. Defaults to None which avoids normalization.
            The key pairs to be used are:
                key 0: run_type {string} -- The options are:
                'shot' - direct shot measurements 
                's005' - 5 second storage
                's020' - 20 second storage
                's100' - 100 second storage
                key 1: parameter {string} -- The options are:
                'N_0'     - counts at time 0
                'y'       - loss rate
                'N_0_err' - associated error
                'y_err'   - associated error
    
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
    """
    
    # instantiate configuration and run type lists
    config_list = ['NOMI', 'JPTI', 'JPSU', 'DISK', 'GD01', 'GD03', 'EPSU']
    run_type_list = ['shot', 's005', 's020', 's100']

    # instantiate the dict
    data_dict = {}

    # the 'all', run_type dicts and the 'all', 'all' dict must be instantiated
    # prior to the main for loops
    for run_type in run_type_list:

        data_dict['all', run_type] = np.empty((0,5), float)

    data_dict['all', 'all'] = np.empty((0,5), float)

    for config in config_list:

        # at the start of each config, initialize the empty 'config', 'all'
        #   array

        data_dict[config, 'all'] = np.empty((0,5), float)

        for run_type in run_type_list:

            arr = load_main(config, run_type, norm_dict)

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

    return data_dict

###############################################################################
###############################################################################

def ucn_yield(data_dict, plotting_flag = False):
    """Analyzes the ucn yield over time.
    
    Arguments:
        data_dict {dict} -- A dictionary of n x 5 data array of the results from 
            loading the run data. The number of rows corresponds to the number 
            of runs loaded. The five columns are:
                0 - the run start time in seconds since the experimental start
                1 - the storage time (0 if direct shot)
                2 - the number of UCN counts
                3 - sqrt(N) error in number of UCN counts
                4 - [day].[run number] of measurement
            The key pairs to be used are:
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
        plotting_flag {boolean, optional} -- Flag to turn plotting on. 
            Defaults to False.
    
    Returns:
        norm_dict {dict} -- dictionary of values of the results from the
            ucn yield analysis. Defaults to None which avoids normalization.
            The key pairs to be used are:
                key 0: run_type {string} -- The options are:
                'shot' - direct shot measurements 
                's005' - 5 second storage
                's020' - 20 second storage
                's100' - 100 second storage
                key 1: parameter {string} -- The options are:
                'N_0'     - counts at time 0
                'y'       - loss rate
                'N_0_err' - associated error
                'y_err'   - associated error
    """
    run_type_list = ['shot', 's005', 's020', 's100']

    if (plotting_flag):

        # for an all-in-one figure
        fig_all, ax_all = plt.subplots()
            
        # for counting loop iterations
        text_y_coord = -0.2

        # for colour consistency in plotting
        ax_c = plt.gca()

    # initializing the dictionary to hold the results
    norm_dict = {}

    for run_type in run_type_list:

        # defining a separate variable for more readable plotting code
        arr = data_dict.copy()['NOMI', run_type]

        # performing a linear fit
        popt, pcov = curve_fit(linear_fit, arr[:,0], arr[:,2], 
                            sigma = arr[:,3], absolute_sigma = True)
        
        # saving the fit results to the dictionary
        norm_dict[run_type, 'N_0']     = popt[0]
        norm_dict[run_type, 'y']       = popt[1]
        norm_dict[run_type, 'N_0_err'] = np.sqrt(np.diag(pcov))[0]
        norm_dict[run_type, 'y_err']   = np.sqrt(np.diag(pcov))[1]

        if (plotting_flag):

            # for colour consistency in plotting
            color = next(ax_c._get_lines.prop_cycler)['color']
            
            # for separate figures
            fig, ax = plt.subplots()
            
            # plotting the data by pre-storage time; separate figures
            ax.errorbar(arr[:,0], arr[:,2], yerr = arr[:,3], fmt = '.',
                        label = run_type, color = color)
            
            # plotting the data by pre-storage time; all on one figure
            ax_all.errorbar(arr[:,0], arr[:,2], yerr = arr[:,3], fmt = '.',
                        label = run_type, color = color)


            ax.plot(arr[:,0], linear_fit(arr[:,0], *popt), color = color);
            ax_all.plot(arr[:,0], linear_fit(arr[:,0], *popt), color = color);
            
            # presentation stuff
            # ax.set_yscale('log')
            ax.set_xlabel('Time Elapsed [s]');
            ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            ax.set_ylabel('UCN Counts');
            ax.legend();
            ax.set_title('Nominal Configuration - Main Detector');

        #     printing the fit results below the figure
        #     text_y_coord = text_y_coord - 0.1
            ax.text(0, text_y_coord, run_type 
                    + r': $N_0 = $%.2e $\pm $ %.2e$, \quad \gamma_{sD_2} = $%.2e $ \pm $ %.2e' % (
                        decimal.Decimal(norm_dict[run_type, 'N_0']), 
                        decimal.Decimal(norm_dict[run_type, 'y']),
                        decimal.Decimal(norm_dict[run_type, 'N_0_err']), 
                        decimal.Decimal(norm_dict[run_type, 'y_err'])),
                    transform=ax.transAxes);

    if (plotting_flag):
        # presentation stuff
        ax_all.set_yscale('log')
        ax_all.set_xlabel('Time Elapsed [s]');
        ax_all.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax_all.set_ylabel('UCN Counts');
        ax_all.legend();
        ax_all.set_title('Nominal Configuration - Main Detector');

    
    return norm_dict

###############################################################################
###############################################################################

def load_p_beam_data(start_time):
    """Loads all of the proton beam data into an array
    
    Arguments:
        start_time {int} -- Unix epoch time stamp of the reference t=0 
            measurement for this data. One should access from the 
            dictionary returned by get_start_times().
    
    Returns:
        p_beam_data (numpy.float64): a (2 x n) array of all of the proton beam 
            monitoring data. 
                - row 0: time elapsed in seconds since the first measurement
                - row 1: beam current in uA
    """
    
    # instantiate array to hold the resulting data, empty and single column 
    # at first, for data to be successively stacked
    p_beam_data = np.empty((0,2), float)
    
    # loop through the files and load the data
    for filename in os.listdir('../data_p_beam_from_edgard'):
        
        # all of the csv file is converted to a list of strings for extracting
        # the time data
        f = open('../data_p_beam_from_edgard/' + filename)
        lines = f.readlines()
        f.close()
        
        # instantiate an array to hold the measurement times
        arr = np.zeros((np.shape(lines)[0], 2))
        
        # loop over every row in the csv file
        for i in range(0, np.shape(lines)[0]):
            
            # convert the measurement time to epoch time
            date_time = lines[i][0:10] + ' ' + lines[i][11:19]
            pattern = '%Y/%m/%d %H:%M:%S'
            measurement_time = int(
                time.mktime(time.strptime(date_time, pattern)))
            
            # save the elapsed time to the times array
            arr[i,0] = measurement_time - start_time
        
        # the current data is loaded into a numpy array
        arr[:,1] = np.loadtxt('../data_p_beam_from_edgard/' + filename, 
                          delimiter = '\t', usecols=(2));

        # removing the 0 values
        for i in range(0,np.shape(arr)[0]):

            if (arr[i,1] == 0):

                arr[i,1] = float('nan')
        
        # append the time and count data to the array
        p_beam_data = np.append(p_beam_data, 
                                arr, 
                                axis = 0)
        
    return p_beam_data

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

# ###############################################################################
# ###############################################################################   

# def load_p_beam_data_3(start_time):
#     """Loads all of the third proton beam data set (p_beam_data_3) into a numpy 
#         array
    
#     Arguments:
#         start_time {int} -- Unix epoch time stamp of the reference t=0 
#             measurement for this data. One should access from the 
#             dictionary returned by get_start_times().
    
#     Returns:
#         p_beam_data (numpy.float64): a (2 x n) array of all of the proton beam 
#             monitoring data. 
#                 - row 0: time elapsed in seconds since the first measurement
#                 - row 1: beam current in uA - monitoring data
#     """
    
#     # instantiate array to hold the resulting data, empty and single column 
#     # at first, for data to be successively stacked
#     p_beam_data = np.empty((0,2), float)
    
#     # loop through the files and load the data
#     for filename in os.listdir('../data_p_beam_3'):
        
#         # all of the csv file is converted to a list of strings for extracting
#         # the time data
#         f = open('../data_p_beam_3/' + filename)
#         lines = f.readlines()
#         f.close()
        
#         # instantiate an array to hold the measurement times
#         arr = np.zeros((np.shape(lines)[0] - 1, 2))
        
#         # loop over every row in the csv file, skipping line 1
#         for i in range(0, np.shape(arr)[0]):
            
#             # convert the measurement time to epoch time
#             date_time = lines[i + 1][0:10] + ' ' + lines[i + 1][11:19]
#             # print(date_time)
#             pattern = '%Y-%m-%d %H:%M:%S'
#             measurement_time = int(
#                 time.mktime(time.strptime(date_time, pattern)))
            
#             # save the elapsed time to the times array
#             arr[i, 0] = measurement_time - start_time

#         # the current data is loaded into a numpy array
#         arr[:,1] = np.loadtxt('../data_p_beam_3/' + filename, 
#                               delimiter = ',', 
#                               skiprows=1, 
#                               usecols=(1));
        
#         # removing the 0 values
#         for i in range(0,np.shape(arr)[0]):

#             if (arr[i,1] == 0):

#                 arr[i,1] = float('nan')


#         # append the time and count data to the array
#         p_beam_data = np.append(p_beam_data, 
#                                 arr, axis = 0)

#     return p_beam_data

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