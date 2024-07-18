import yaml, glob
import xarray as xr
import numpy as np



def _format_time (dataset, date_start, date_end, freq, leap_years=True, time_label='time'):
    """ 
    Function that takes a dataset with a time coordinate and changes the time coordinate format to a datetime64[ns] format.

    Args:
        dataset:            the dataset with the time coordinate
        date_start:         the date of the first measurement. 'YYYY-MM-DD' for daily, 'YYYY-MM' for monthly, etc.
        date_end:           the date of the last measurement. 'YYYY-MM-DD' for daily, 'YYYY-MM' for monthly, etc.
        freq:               frequency of measurements. 'D' for daily, 'M' for monthly, 'Y' for yearly.
        leap_years:         if False, then all Februaries have 28 days.      
        time_label:         str label of the 'time' coordinate in the dataset
        
    Returns:
        dataset with the time coordinate reformatted as datetime64[ns].
    """

    if type(freq) == type(None):
        raise TypeError("Not a valid value of 'freq' for formatting of time. Please set 'format_time' to False in the config, or specify a valid frequency.")

    date_start = np.datetime64(date_start)
    date_end = np.datetime64(date_end)
    date_range = np.arange(date_start, date_end + np.timedelta64(1,freq), np.timedelta64(1,freq)).astype('datetime64[ns]')
    
    # Deal with leap years, if necessary.
    if freq == 'D':
        # Remove ALL February 29ths from date_range
        if leap_years == False:
            date_range = date_range[np.array([not (date[5:10] == '02-29') for date in date_range.astype(str)])]
        # Else, remove all February 29ths from date_range UNLESS it's a leap year.
        else:
            date_range = date_range[np.array([not (date[5:10] == '02-29' and int(date[0:4])%4 != 0) for date in date_range.astype(str)])]

    if not len(dataset[time_label]) == len(date_range):
        raise ValueError(f"Number of points in specified range of time (N={len(date_range)}) is incompatible with the number of data points (N={len(dataset[time_label])}) and the given frequency.")
    dataset[time_label] = date_range
    return dataset



def get_dataset(dataset_label, config_path=None):
    """
    Loads and returns a .nc dataset into memory using the information in config.yaml

    Args:
        dataset_label (str):        The label given to your dataset in config.yaml
        config_path (str):          The path to config.yaml. If None, then assumed to be in pwd.
        
    Returns:
        An xarray dataset, or a list of xarray datasets, depending on how many files there are for the dataset.
    """

    # Open config
    if type(config_path) == type(None):
        config_path = 'config.yaml'
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    config = config['Datasets'][dataset_label]

    # Open dataset from config-specified path. If config specifies multiple paths, then the data for each path is opened and appended to a list.
    paths = sorted(glob.glob(config['path']))
    dataset = xr.open_dataset(paths[0]) if len(paths) == 1 else [xr.open_dataset(path) for path in paths]

    # Format time coordinate, if 'format_time' in the config is true for this dataset.
    if 'format_time' in config:
        if config['format_time'] == True:
            time_settings = config['time']
            if 'leap_years' in time_settings:
                leap_years = time_settings['leap_years']
            else:
                leap_years = True
            time_label = time_settings['label']
            if len(paths) == 1:
                dataset = _format_time (dataset, 
                                        date_start = time_settings['date_begin'], 
                                        date_end = time_settings['date_end'], 
                                        freq = time_settings['freq'], 
                                        leap_years = leap_years, 
                                        time_label = time_label)
            else:
                dataset = [_format_time (data, 
                                        date_start = time_settings['date_begin'], 
                                        date_end = time_settings['date_end'], 
                                        freq = time_settings['freq'], 
                                        leap_years = leap_years, 
                                        time_label = time_label) for data in dataset]
    return dataset
