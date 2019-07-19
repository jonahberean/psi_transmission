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
    return N_0 + y * t

###############################################################################
###############################################################################


def get_start_times():
    """Retrieves the UNIX epoch time stamp of various relevant start times
    
    Returns:
        dict: the UNIX epoch time stamps in a dictionary. Key options are:
            'edgard_p_beam' - first proton beam current measurement (edgard)
            'bernhard_p_beam' - first proton beam current measurement 
                (bernhard)
    """
    
    # initialize the dictionary
    start_time_dict = {}
    
    ## 'edgard_p_beam': Proton beam current
    
    # read the 8th of December data as a list of strings
    f = open('../data_p_beam_from_edgard/20171208.csv')
    lines = f.readlines()
    f.close()
    
    # convert the measurement time to epoch time
    date_time = lines[0][0:10] + ' ' + lines[0][11:19]
    pattern = '%Y/%m/%d %H:%M:%S'
    measurement_time = int(time.mktime(time.strptime(date_time, pattern)))
    start_time_dict['edgard_p_beam'] = measurement_time

    ## 'bern_p_beam': Proton beam current from Bernhard

    # read the 8th of December data as a list of strings
    f = open('../data_p_beam/2cbb80_analog_20171208_0007.csv')
    lines = f.readlines()
    f.close()
    
    # convert the measurement time to epoch time
    date_time = lines[1][0:10] + ' ' + lines[1][11:19]
    pattern = '%d.%m.%Y %H:%M:%S'
    measurement_time = int(time.mktime(time.strptime(date_time, pattern)))
    start_time_dict['bernhard_p_beam'] = measurement_time
    
    return start_time_dict

###############################################################################
###############################################################################

def load_data(config, run_type, start_time, norm_dict = None):
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
                's000' - direct shot measurements 
                's005' - 5 second storage
                's020' - 20 second storage
                's100' - 100 second storage
        start_time {int} -- Unix epoch time stamp of the reference t=0 
            measurement for this data. One should access from the 
            dictionary returned by get_start_times().
        norm_dict {dict} -- dictionary of values of the results from the
            ucn yield analysis. Defaults to None which avoids normalization.
            The key pairs to be used are:
                key 0: run_type {string} -- The options are:
                's000' - direct shot measurements 
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
        
    # instantiate a new numpy array 
    arr = np.empty((0,5), float)
            
    for filename in os.listdir('../data_main/sorted'):

        # Only the files matching our desired configuration and run 
        # type are selected. The '.tof' condition is just so we 
        # don't perform the analysis twice per run (since it would
        # otherwise match to the .tof and the .txt files)
        if ((config in filename) and (run_type in filename) and 
        ('.tof' in filename)):

            # grab from the text file associated with the run
            f = open( '../data_main/sorted/' + filename[0:22] 
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
            storage_time = int(run_type[1:4])

            # The data is retrieved from the .tof file
            count_data = np.loadtxt('../data_main/sorted/' + 
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

                counts = np.sum(count_data[150:1000])

            # specific data cut for run 66 on the 8th
            elif ((filename[2:3] == '8') and 
                  (filename[10:12] == '66')):

                counts = np.sum(count_data[150:1500])

            # specific data cut for run 88 on the 8th
            elif ((filename[2:3] == '8') and 
                  (filename[10:12] == '88')):

                counts = np.sum(count_data[150:2500])

            # if it is a shot run we take all the counts
            # !!! But does this mean that we're counting that initial proton 
            # irradiation background in here? What is that background from?
            # should it be chracterized and removed? Does it fluctuate over
            # time? If so could it be tied to the current measurements and then
            # corrected out?
            elif (run_type == 's000'):
            
                counts = np.sum(count_data)
            
            # otherwise cut the data normally for a pre-storage run
            # this cuts out the initial background appearing from
            # irradiation
            else:

                counts = np.sum(count_data[150:-1])

            # normalize the data depending on the normalize_flag
            if (norm_dict != None):

                # !!! the error associated with the normalization 
                # routine has not been implemented

                # the fits to the nominal configuration data provide
                # a benchmark for percentage loss of absolute counts
                # depending on the time having elapsed since the start
                # of the experiment. 

                # interpolate on the nominal fit to get a counts value
                # for this run time
                interp_counts = linear_fit(run_time, 
                                           norm_dict[run_type, 'N_0'], 
                                           norm_dict[run_type, 'y'])

                # the interpretation here is that whatever count total 
                # the current run has, it's total counts are actually 
                # this norm_factor of what they would have been without
                # the degradation of the sD2 surface.
                norm_factor = interp_counts / norm_dict[run_type,
                                                        'N_0']

                # this brings the counts back to what they would have 
                # been
                counts = counts / norm_factor

            # saving the [day].[run number] can be useful for debugging
            day_run_no = int(filename[1:3]) + (0.001
                                               * int(filename[9:12]))

            # append the loaded data to the existing array
            arr = np.append(arr, [[run_time,
                                   storage_time,
                                   counts,
                                   np.sqrt(counts),
                                   day_run_no]],
                            axis=0)

    # return the completed array, sorted along the time axis
    return arr[arr[:,0].argsort()]

###############################################################################
###############################################################################

def load_all_data(start_time, norm_dict = None):
    """A function to load data and sum counts for all the run data available
    
    Arguments:
        start_time {int} -- Unix epoch time stamp of the reference t=0 
            measurement for this data. One should access from the 
            dictionary returned by get_start_times().
        norm_dict {dict} -- dictionary of values of the results from the
            ucn yield analysis. Defaults to None which avoids normalization.
            The key pairs to be used are:
                key 0: run_type {string} -- The options are:
                's000' - direct shot measurements 
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
                's000' - direct shot measurements 
                's005' - 5 second storage
                's020' - 20 second storage
                's100' - 100 second storage
                'all'  - all of the above
    """
    
    # instantiate configuration and run type lists
    config_list = ['NOMI', 'JPTI', 'JPSU', 'DISK', 'GD01', 'GD03', 'EPSU']
    run_type_list = ['s000', 's005', 's020', 's100']

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

            arr = load_data(config, run_type, start_time, norm_dict)

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
                's000' - direct shot measurements 
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
                's000' - direct shot measurements 
                's005' - 5 second storage
                's020' - 20 second storage
                's100' - 100 second storage
                key 1: parameter {string} -- The options are:
                'N_0'     - counts at time 0
                'y'       - loss rate
                'N_0_err' - associated error
                'y_err'   - associated error
    """
    # list of run types
    run_type_list = ['s000', 's005', 's020', 's100']
    
    # initializing the dictionary to hold the results
    norm_dict = {}
    
    # for colour consistency in plotting
    ax = plt.gca()

    # for counting loop iterations
    text_y_coord = -0.2

    for run_type in run_type_list:
        
        # defining a separate variable for more readable plotting code
        arr = data_dict['NOMI', run_type]

        # for colour consistency in plotting
        color = next(ax._get_lines.prop_cycler)['color']

        # plotting the data by pre-storage time
        plt.errorbar(arr[:,0], arr[:,2], yerr = arr[:,3], fmt = '.',
                     label = run_type, color = color)

        # performing a linear fit
        popt, pcov = curve_fit(linear_fit, arr[:,0], arr[:,2], 
                               sigma = arr[:,3], absolute_sigma = True)
        plt.plot(arr[:,0], linear_fit(arr[:,0], *popt), color = color);

        # saving the fit results to the dictionary
        norm_dict[run_type, 'N_0']     = popt[0]
        norm_dict[run_type, 'y']       = popt[1]
        norm_dict[run_type, 'N_0_err'] = np.sqrt(np.diag(pcov))[0]
        norm_dict[run_type, 'y_err']   = np.sqrt(np.diag(pcov))[1]

        # printing the fit results below the figure
        text_y_coord = text_y_coord - 0.1
        plt.text(0, text_y_coord, run_type 
                 + r': $N_0 = $%.2e $\pm $ %.2e$, \quad \gamma_{sD_2} = $%.2e $ \pm $ %.2e' % (
                     decimal.Decimal(norm_dict[run_type, 'N_0']), 
                     decimal.Decimal(norm_dict[run_type, 'y']),
                     decimal.Decimal(norm_dict[run_type, 'N_0_err']), 
                     decimal.Decimal(norm_dict[run_type, 'y_err'])),
                 transform=ax.transAxes);

    # presentation stuff
    plt.yscale('log')
    plt.xlabel('Time Elapsed [s]');
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ylabel('UCN Counts');
    plt.legend();
    plt.title('Nominal Configuration - Main Detector');
    
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

def load_bern_p_beam_data(start_time):
    """Loads all of the proton beam data (Bernhard's) into an array
    
    Arguments:
        start_time {int} -- Unix epoch time stamp of the reference t=0 
            measurement for this data. One should access from the 
            dictionary returned by get_start_times().
    
    Returns:
        p_beam_data (numpy.float64): a (2 x n) array of all of the proton beam 
            monitoring data. 
                - row 0: time elapsed in seconds since the first measurement
                - row 1: beam current in uA - monitoring data
                - row 3: beam current in uA - timing data (not to be trusted
                    for absolute value)
    """
    
    # instantiate array to hold the resulting data, empty and single column 
    # at first, for data to be successively stacked
    p_beam_data = np.empty((0,3), float)
    
    # loop through the files and load the data
    for filename in os.listdir('../data_p_beam'):
        
        # all of the csv file is converted to a list of strings for extracting
        # the time data
        f = open('../data_p_beam/' + filename)
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
        arr[:,1:3] = np.loadtxt('../data_p_beam/' + filename, 
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

def load_p_beam_data_3(start_time):
    """Loads all of the third proton beam data set (p_beam_data_3) into a numpy 
        array
    
    Arguments:
        start_time {int} -- Unix epoch time stamp of the reference t=0 
            measurement for this data. One should access from the 
            dictionary returned by get_start_times().
    
    Returns:
        p_beam_data (numpy.float64): a (2 x n) array of all of the proton beam 
            monitoring data. 
                - row 0: time elapsed in seconds since the first measurement
                - row 1: beam current in uA - monitoring data
    """
    
    # instantiate array to hold the resulting data, empty and single column 
    # at first, for data to be successively stacked
    p_beam_data = np.empty((0,2), float)
    
    # loop through the files and load the data
    for filename in os.listdir('../data_p_beam_3'):
        
        # all of the csv file is converted to a list of strings for extracting
        # the time data
        f = open('../data_p_beam_3/' + filename)
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
        arr[:,1] = np.loadtxt('../data_p_beam_3/' + filename, 
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

def find_coincidences(p_beam_data, main_data_dict, plotting_flag = False):
    """Searches for coincident runs and measurements, from the main detector
        data and the proton beam current data, respectively.
    
    Arguments:
        p_beam_data (numpy.float64) -- a (2 x n) array of all of the proton beam 
            monitoring data. 
                - row 0: time elapsed in seconds since the first measurement
                - row 1: beam current in uA - monitoring data
        main_data_dict (dict) -- see docstring of load_all_data()
        plotting_flag (boolean, optional) -- flag to enable plotting. Defaults
            to False.
    
    Returns:
        (dict) -- a dictionary of arrays for coincident proton beam
            measurements, labelled by run type. The key options are:
            key 0: run_type {string} -- The options are:
                's000' - direct shot measurements 
                's005' - 5 second storage
                's020' - 20 second storage
                's100' - 100 second storage
    """

    # instantiate a new dictionary
    reduced_dict = {}

    # iterate over run types
    # run_type_list = ['s000', 's005', 's020', 's100']
    run_type_list = ['s000']
    for run_type in run_type_list:

        n_arr = main_data_dict['all', run_type]
        p_arr = p_beam_data

        # instantiate a new array to hold the reduced data set
        # red_arr = np.empty((0,2), float)

        # for i in range(0, np.shape(n_arr)[0]):
        for i in range(0, 1):

            red_arr = np.where((np.abs(n_arr[i,0] - p_arr[:,0]) < 20), 
                                p_arr[:,1], 
                                np.zeros(np.shape(p_arr)[0]))

            red_arr = red_arr[np.nonzero(red_arr)]  

            print(red_arr)
            
            red_arr = np.where(np.isnan(red_arr),
                                np.zeros(np.shape(red_arr)[0]),
                                red_arr)

            red_arr = red_arr[np.nonzero(red_arr)]    

            print(red_arr)                
            # condition = np.abs(p_arr[:,0] - n_arr[i,0]) < 20
            # condition


        # # iterate over every run start time, and every proton beam measurement time
        # # do this in a nested format
        # for i in range(0, np.shape(n_arr)[0]):
        
        #     for j in range(0, np.shape(p_arr)[0]):
        #         # check if non-zero
        #         if (p_arr[j,0] != float('nan')):
        #             # calculate the absolute time difference between the two times
        #             abs_diff = abs(n_arr[i,0] - p_arr[j,0])
                    
        #             # calculate the non-absolute time difference between the two times
        #             diff = n_arr[i,0] - p_arr[j,0]
                    
        #             # if the measurement time was within 9 seconds of the run time
        #             if (abs_diff < 9):

        #                 red_arr = np.append(red_arr, [p_arr[j,:]], axis = 0)
        
        # # store the completed array in the dictionary
        # reduced_dict[run_type] = red_arr

        # if (plotting_flag):

        #         # instantiate the subplots 
        #         fig, ax1 = plt.subplots()

        #         # plot the proton beam current data
        #         ax1.scatter(red_arr[:,0], red_arr[:,1], s=1, color = 'r')

        #         # presentation stuff
        #         ax1.set_xlabel('Time Elapsed (s)')
        #         ax1.set_ylabel(r'Proton Beam Current [$\mu$A]', 
        #             color = 'r')
        #         ax1.tick_params(axis='y')
        #         ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

        #         # instantiate a second axes that shares the same x-axis
        #         ax2 = ax1.twinx()  
                
        #         # plot the neutron count data for runs of this pre-storage
        #         # time
        #         ax2.errorbar(n_arr[:,0], n_arr[:,2], yerr = n_arr[:,3], 
        #                         fmt = '.', color = 'b')

        #         # presentation stuff
        #         ax2.set_ylabel('UCN Counts')  
        #         ax2.set_yscale('log')
        #         fig.tight_layout()  

    # return reduced_dict
    return red_arr


###############################################################################