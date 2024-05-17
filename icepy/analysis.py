### ---------------------------------------------
### Author: Robert Payne
### Contact: gpayne1654@uvic.ca
###
### This package makes extensive use of array indexing techniques, and the libraries numpy, pandas, xarray, and matplotlib. 
### It's recommended you familiarize yourself with these before working with this library.
### ---------------------------------------------

import numpy as np
import xarray as xr
import pandas as pd
import cftime

def _cftime_to_datetime64(xr_data):
    """
    Convert 'time' coordinate of an xarray Dataset or DataArray from cftime objects to datetime64[ns] objects.
    Not callable external to the module.
    
    Args:
        xr_data (xarray.Dataset or xarray.DataArray): The xarray Dataset or DataArray containing a 'time' coordinate.
        
    Returns:
        the input dataset, with the 'time' coordinate as datetime64[ns] object
    """
    # Extract the 'time' coordinate from the xarray object
    time_coord = xr_data['time']
    
    # Check if 'time' coordinate is already in pandas datetime format
    if isinstance(time_coord.values[0], pd.Timestamp):
        return time_coord
    
    # Check if 'time' coordinate is in cftime format
    elif isinstance(time_coord.values[0], cftime.datetime):
        # Convert cftime dates to pandas datetimes
        pd_datetimes = []
        for date in time_coord.values:
            if isinstance(date, cftime.datetime):
                pd_datetimes.append(pd.Timestamp(date.year, date.month, date.day, date.hour, date.minute, date.second, date.microsecond))
            else:
                raise TypeError("Unsupported time coordinate format. Supported formats are pandas datetimes and cftime datetimes.")
        xr_data['time'] = pd.DatetimeIndex(pd_datetimes)
        return xr_data
    
    else:
        raise TypeError("Unsupported time coordinate format. Supported formats are pandas datetimes and cftime datetimes.")



def for_month (dataset, month, drop=True):
    """
    Function that takes in a dataset with a time coordinate and returns a dataset with only the data for the specified month.

    Args:
        dataset:                the dataset with the time coordinate
        month (int or list):    the month (or list of months) for which you want the data
        drop (bool):            whether to drop the data that aren't for the given month. If False, then these data are set to nans instead. Default True.

    Returns:
        dataset with the data for the corresponding month(s).
    """

    if type(month) == int:
        return dataset.where(dataset['time.month'] == month, drop=drop)
    else:
        return dataset.where(dataset['time.month'].isin(month), drop=drop)



def format_time_coord (dataset, date_start, date_end, freq):
    """ 
    Function that takes a dataset with a time coordinate and changes the time coordinate format to a datetime64[ns] format.

    Args:
        dataset:            the dataset with the time coordinate
        date_start:         the date of the first measurement. 'YYYY-MM-DD' for daily, 'YYYY-MM' for monthly, etc.
        date_end:           the date of the last measurement. 'YYYY-MM-DD' for daily, 'YYYY-MM' for monthly, etc.
        freq:               frequency of measurements. 'D' for daily, 'M' for monthly, 'Y' for yearly.

    Returns:
        dataset with the time coordinate reformatted as datetime64[ns].
    """

    date_start = np.datetime64(date_start)
    date_end = np.datetime64(date_end)
    date_range = np.arange(date_start, date_end + np.timedelta64(1,freq), np.timedelta64(1,freq))
    if not len(dataset['time']) == len(date_range):
        raise ValueError("Specified range of time is incompatible with the number of data points and the given frequency.")
    dataset['time'] = date_range
    return dataset



def sic_to_sie (sic_dataset, grid_area_dataset, lat_bounds=None, lon_bounds=None, lat_label='lat', lon_label='lon', sic_label='SICN', flip_meridional=False, flip_zonal=False, ensemble=None, mfactor=1e-12, as_dataframe=False):
    """
    Function that takes in a netcdf dataframe of gridded sea ice concentration (SIC) and outputs the corresponding sea ice extent (in units of 10^6 km^2) 
    at each time step, conventionally defined as the integrated area of all grid cells having SIC > 0.15

    Args:
        sic_dataset:                an xarray or netcdf dataset containing a SIC variable as a function of coordinates latitude, longitude, and time.
        grid_area_dataset:          an xarray or netcdf dataset containing grid cell area as a function of coordinates latitude and longitude.
                                    Can be computed through the climate data operators command 'cdo gridarea'.
        lat_bounds (tuple):         the lower and upper bound of latitudes in which SIC will be considered. Outside of the bounds, SIC is set to nan.
        lon_bounds (tuple):         similar to lat_bounds, except for longitude. Check whether your dataset inherently uses (0,360) or (-180,180).
        lat_label (string):         label of the latitude coordinate in the dataset (not the grid file). Often 'lat' or 'latitude'. 
        lon_label (string):         label of the longitude coordinate in the dataset (not the grid file). Often 'lon' or 'longitude'.
        sic_label (string):         label of the SIC data variable in the dataset. Often 'SICN' or 'siconc'.
        flip_meridional (bool):     if True, then sic INSIDE the LATITUDE bounds are set to nans, rather than outside. Defaults to False.
        flip_zonal (bool):          if True, then sic INSIDE the LONGITUDE bounds are set to nans, rather than outside. Defaults to False.
        ensemble (None, str, int):  selects the given ensemble and drops the rest. if 'ave' or 'mean', then averages over ensembles instead.
        mfactor (float):            multiplicative factor to multiple the final result by. Use 1e-12 if you want the outputted SIE to be in units of 10^6 sq km. If 1, then output is in units of sq m.
        as_dataframe (bool):        if 'get' and this are true, then this function will return the SIE dataset as a pandas dataframe rather than an xarray dataset. Defaults to False.

    Returns:
        (xarray dataset, optional)  the dataset of calculated SIE.
    """

    grid_area = grid_area_dataset
    SIC = sic_dataset
    if grid_area.dims == {'lon':360, 'lat':180}:
        grid_area = grid_area.rename({'lat':lat_label,'lon':lon_label})

    # apply latitude and longitude bounds
    if not lat_bounds == None:
        lat_lo, lat_up = lat_bounds
        if lat_lo > lat_up:
            raise ArithmeticError("Lower latitude bound cannot be greater than upper latitude bound.")
        if flip_meridional:
            SIC = SIC.where((SIC[lat_label] <= lat_lo) | (SIC[lat_label] >= lat_up)) # OR
        else:
            SIC = SIC.where((SIC[lat_label] >= lat_lo) & (SIC[lat_label] <= lat_up)) # AND
    if not lon_bounds == None:
        lon_lo, lon_up = lon_bounds
        if lon_lo > lon_up:
            raise ArithmeticError("Lower longitude bound cannot be greater than upper longitude bound.")
        if flip_zonal:
            SIC = SIC.where((SIC[lon_label] <= lon_lo) | (SIC[lon_label] >= lon_up)) # OR
        else:
            SIC = SIC.where((SIC[lon_label] >= lon_lo) & (SIC[lon_label] <= lon_up)) # AND

    # calculate SIE
    SIE = grid_area.expand_dims(time=SIC['time']).where(SIC[sic_label] >= 0.15).sum(dim=(lat_label,lon_label))
    SIE *= mfactor
    # if method == 'old':  *****DONT USE THIS, WILL BE REMOVED AT SOME POINT, IT'S LESS ACCURATE AND MORE COMPUTATIONALLY EXPENSIVE
    #     SIE = grid_area.expand_dims(time=SIC['time']).where(SIC[sic_label] >= 0.15)
    #     SIE *= (111120**2)*np.abs(np.cos(SIE.lat*np.pi/180))
    #     SIE = SIE.sum(dim=(lat_label,lon_label))
    #     SIE *= mfactor
    SIE = SIE.rename({'cell_area': 'SIE'})

    # average/select ensembles if specified
    if not ensemble == None:
        if ensemble == 'average' or ensemble == 'ave' or ensemble == 'mean':
            SIE = SIE.groupby('ensemble', squeeze=False).mean()
        else:
            SIE = SIE.where(SIE['ensemble'] == ensemble, drop=True)
        SIE['SIE'] = SIE['SIE'][list(SIE.dims)[::-1].index('ensemble')]  # removes ensemble as a dimension (xr.drop_dims obliterates SIE, for some reason)

    # match time coordinate format with that of input and then return the dataset
    SIE['time'] = SIC['time']
    if as_dataframe:
        return SIE.to_dataframe()
    else:
        return SIE



def get_climatology (dataset, var, ref_period=None):
    """
    Function that takes in a netcdf dataset and calculates the climatology for a user-specified variable and reference period.

    Args:
        dataset (xarray dataset):   the dataset containing a time coordinate and the variable for which the climatology is to be calculated for.
        var (string):               the label of the variable for the climatology (e.g., "SIE")
        ref_period (tuple):         a tuple of two years defining the timespan (inclusive) over which the climatology will be calculated. If 'None', then entire time span is used.

    Returns:
        (1x13 array):   The climatology. The first entry is the mean over all months, whereas the remaining twelve correspond to the mean for each specific month, starting in January.
    """

    # include only the data within the reference period range
    if not ref_period == None:
        start_year, end_year = ref_period
        dataset = dataset.where((dataset['time.year'] >= start_year) & (dataset['time.year'] <= end_year))

    # calculate the climatology for each month
    for i in range(13):
        if i == 0:
            mean_over_all_months = [dataset.mean(dim='time')[var].values]
        else:
            means_for_each_month = dataset.groupby('time.month').mean(dim='time')[var].values
    climatology = np.concatenate((mean_over_all_months, means_for_each_month))

    return climatology



def get_anomalies (dataset, var, ref_period=None):
    """
    Function that takes in a netcdf dataset and calculates the climatology (based on a mean) for a user-specified variable and reference period.

    Args:
        dataset (xarray dataset):   the dataset containing a time coordinate and the variable from which anomalies will be calculated.
        var (string):               the label of the variable to calculate anomalies (e.g., "SIE")
        ref_period (tuple):         a tuple of two years defining the timespan (inclusive) over which the climatology will be calculated. If 'None', then entire time span is used.

    Returns:
        (dataset, same dims as input): dataset as anomalies, with respect to the climatology calculated over the reference period provided.
    """

    # climatology for all months
    climatology = get_climatology(dataset=dataset, var=var, ref_period=ref_period)[1:]
    
    # subtract climatology from dataset to produce anomalies
    dataset[var] = dataset[var].groupby('time.year').map(lambda x: x - climatology[x['time.month'].values[0]-1 : x['time.month'].values[-1]])
    return dataset



def get_trend (dataset):
    """
    Function that computes a least-squares linear fit for the given data variable, assuming the data are equidistant in time.
    Currently this function is only compatible with a 1-dimensional time series. (i.e., the only coordinate is time.)

    Args:
        dataset (xarray dataset): the dataset containing a time coordinate and the variable to be detrended.

    Returns:
        (1D list): the trend line values.
    """

    N = len(dataset['time'])
    slope, intercept = np.polyfit(np.arange(0,N,1), dataset, deg=1)
    return slope*np.arange(0,N,1) + intercept



def remove_trend (dataset, var, method='mean', ref_period=(0,9999)):
    """
    Function that removes linear trends using Mitchell Bushuk's method.
    In short, the trend that's removed from the ith data point in the time series is calculated using only the first i data points (i, i-1, ..., 2, 1).
    For i = 1,2,3, the mean is removed instead.
    This is done separately for every month (i.e., the trend calculated for a november data point only uses other data pts in november for its calculation).
    This function is only compatible with a 1-dimensional time series. (i.e., the only coordinate is time.) (unless method='mean', in which case it should always work.)

    Args:
        dataset (xarray dataset):   the dataset containing a time coordinate and the variable to be detrended.
        var (string):               the label of the variable to calculate anomalies (e.g., "SIE")
        method (string):            if 'Mitch', then use Mitch's detrending method. If 'Mean', then remove mean. If 'linear', just a standard linear detrending.
        ref_period (tuple):          use if you want the trend to be calculated only using data from a range of years, e.g., between the years 2000 and 2020.

    Returns:
        (dataset, same dims as input): dataset with the trend removed.
    """

    dataset = dataset.copy(deep=True) # deep must be set to True

    for m in np.arange(1,13,1): # loop over each month

        # get subset of data for that month
        month_indices = (dataset['time.month'] == m) & (dataset['time.year'] >= ref_period[0]) & (dataset['time.year'] <= ref_period[1])
        subset = dataset.where(month_indices,drop=True)
        subset_detrended = subset.copy(deep=True) # deep must be set to True

        # Mitch's detrending method
        if method == 'Mitch' or method == 'mitch': 

            for i in range(len(subset[var])): # loop over each point for that month and apply the detrending procedure

                # for i = 0,1,2, remove the mean
                if i <= 2:
                    
                    mean_for_month = subset[var][:i+1].mean()
                    subset_detrended[var][i] -= mean_for_month

                # for i > 2, remove the corresponding linear trend
                else:

                    trend_for_month = get_trend(subset[var][:i+1])
                    subset_detrended[var][i] -= trend_for_month[i]

            dataset[var].loc[month_indices] = subset_detrended[var]

        # Just remove the mean
        elif method == 'mean':
            subset_detrended[var] = subset[var] - subset[var].mean(dim='time')
            dataset[var].loc[month_indices] = subset_detrended[var]

        # "Standard" detrending method
        elif method == 'linear' or method == 'lin':
            subset_detrended[var] = subset[var] - get_trend(subset[var])
            dataset[var].loc[month_indices] = subset_detrended[var]

        else:
            raise ValueError("Not a valid detrending method. Choose one of ['mean', 'linear', 'Mitch']")

    return dataset


