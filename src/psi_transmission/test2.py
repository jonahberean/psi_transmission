###############################################################################

def load_data(config, run_type, start_point, norm_dict = None):
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
        start_point {int} -- Unix epoch time stamp of the reference t=0 
            measurement for this data
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

            # cut the data normally
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