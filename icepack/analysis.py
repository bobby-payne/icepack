### ---------------------------------------------
### Author: Robert Payne
### Contact: gpayne1654@uvic.ca
###
### This package makes extensive use of array indexing techniques, and the libraries numpy, xarray, and matplotlib. 
### It's recommended you familiarize yourself with these before working with this library.
### ---------------------------------------------

import numpy as np
import xarray as xr
import pandas as pd
import cftime



def _apply_bounds(dataset, lat_bounds, lon_bounds, lat_label, lon_label, flip_zonal=False, flip_meridional=False):
    """
    Apply latitude and longitude bounds to a dataset, and return the dataset with all variables outside the bounds set to nans.

    Args:
        dataset:                    a dataset containing a variable as a function of coordinates latitude, longitude to which to apply the bounds.
        lat_bounds (tuple):         the lower and upper bound of latitudes in which the variable will be considered. Outside of the bounds, the variable is set to nan.
        lon_bounds (tuple):         similar to lat_bounds, except for longitude. Check whether your dataset inherently uses (0,360) or (-180,180).
        lat_label (string):         label of the latitude coordinate in the dataset (not the grid file). Often 'lat' or 'latitude'. 
        lon_label (string):         label of the longitude coordinate in the dataset (not the grid file). Often 'lon' or 'longitude'.
        flip_meridional (bool):     if True, then sic INSIDE the LATITUDE bounds are set to nans, rather than outside. Defaults to False.
        flip_zonal (bool):          if True, then sic INSIDE the LONGITUDE bounds are set to nans, rather than outside. Defaults to False.
    """

    if not lat_bounds == None:
        lat_lo, lat_up = lat_bounds
        if lat_lo > lat_up:
            raise ArithmeticError("Lower latitude bound cannot be greater than upper latitude bound.")
        if flip_meridional:
            dataset = dataset.where((dataset[lat_label] <= lat_lo) | (dataset[lat_label] >= lat_up)) # OR
        else:
            dataset = dataset.where((dataset[lat_label] >= lat_lo) & (dataset[lat_label] <= lat_up)) # AND
    if not lon_bounds == None:
        lon_lo, lon_up = lon_bounds
        if lon_lo > lon_up:
            raise ArithmeticError("Lower longitude bound cannot be greater than upper longitude bound.")
        if flip_zonal:
            dataset = dataset.where((dataset[lon_label] <= lon_lo) | (dataset[lon_label] >= lon_up)) # OR
        else:
            dataset = dataset.where((dataset[lon_label] >= lon_lo) & (dataset[lon_label] <= lon_up)) # AND

    return dataset



def select_date(dataset,year=None,month=None,day=None,drop=True):
    """
    Selects the data of a particular xarray dataset for a particular year, month, and day.

    Args:
        dataset:                the dataset with the time coordinate
        year (int or None):     the year for which you want the data
        month (int or None):    the month for which you want the data
        day (int or None):      the day for which you want the data
        drop (bool):            whether or not to drop the data outside the given date. Default True.

    Returns:
        dataset with the data for the given year, month, and/or day.
    """

    data = dataset
    if not year == None:
        try:
            data = data.where(data['time.year']==year,drop=drop)
        except:
            data = data.where(data['year']==year,drop=drop)
    if not month == None:
        try:
            data = data.where(data['time.month']==month,drop=drop)
        except:
            data = data.where(data['month']==month,drop=drop)
    if not day == None:
        try:
            data = data.where(data['time.day']==day,drop=drop)
        except:
            data = data.where(data['day']==day,drop=drop)

    return data
    


def format_time_coord (dataset, date_start, date_end, freq, leap_years=True, time_label='time'):
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



def get_iceextent (sic_dataset, grid_area_dataset=None, lat_bounds=None, lon_bounds=None, lat_label='lat', lon_label='lon', sic_label='SICN', flip_meridional=False, flip_zonal=False, ensemble=None, ensemble_label='ensemble', sic_factor=1, mfactor=1e-12):
    """
    Function that takes in a netcdf dataframe of gridded sea ice concentration (SIC) and outputs the corresponding sea ice extent (in units of 10^6 km^2, unless specified by user w/ mfactor) 
    at each time step, conventionally defined as the integrated area of all grid cells having SIC > 0.15

    Args:
        sic_dataset:                an xarray or netcdf dataset containing a SIC variable as a function of coordinates latitude, longitude, and time.
        grid_area_dataset:          an xarray or netcdf dataset containing grid cell area as a function of coordinates latitude and longitude.
                                    Can be computed through the climate data operators command 'cdo gridarea'.
                                    If no grid is provided (i.e., left as None), then SIE will be calculated with a less accurate method specific to 1x1 grids.
        lat_bounds (tuple):         the lower and upper bound of latitudes in which SIC will be considered. Outside of the bounds, SIC is set to nan.
        lon_bounds (tuple):         similar to lat_bounds, except for longitude. Check whether your dataset inherently uses (0,360) or (-180,180).
        lat_label (string):         label of the latitude coordinate in the dataset (not the grid file). Often 'lat' or 'latitude'. 
        lon_label (string):         label of the longitude coordinate in the dataset (not the grid file). Often 'lon' or 'longitude'.
        sic_label (string):         label of the SIC data variable in the dataset. Often 'SICN' or 'siconc'.
        flip_meridional (bool):     if True, then sic INSIDE the LATITUDE bounds are set to nans, rather than outside. Defaults to False.
        flip_zonal (bool):          if True, then sic INSIDE the LONGITUDE bounds are set to nans, rather than outside. Defaults to False.
        ensemble (None, str, int):  selects the given ensemble and drops the rest. if 'ave' or 'mean', then averages over ensembles instead.
        ensemble_label (str):       label of the ensemble coordinate in the dataset. Often 'ensemble' or 'realization'.
        sic_factor (float):         multiplies all sic values by this factor. useful for converting from percents. default 1.
        mfactor (float):            multiplicative factor to multiple the final result by. Use 1e-12 if you want the outputted SIE to be in units of 10^6 sq km. If 1, then output is in units of sq m.

    Returns:
        (xarray dataset, optional)  the dataset of calculated SIE.
    """

    grid_area = grid_area_dataset
    SIC = sic_dataset
    if not type(grid_area_dataset) == type(None):
        if grid_area.dims == {'lon':360, 'lat':180}:
            grid_area = grid_area.rename({'lat':lat_label,'lon':lon_label})

    # apply latitude and longitude bounds
    SIC = _apply_bounds(SIC, lat_bounds, lon_bounds, lat_label, lon_label, flip_zonal, flip_meridional)

    # calculate SIE
    SIC[sic_label] *= sic_factor
    SIC = SIC.where((SIC[sic_label] <= 1.)&(SIC[sic_label] >= 0.),drop=True)
    if type(grid_area_dataset) == type(None):
        SIE = SIC.where(SIC[sic_label] >= 0.15,drop=True)
        SIE *= (111120**2)*np.abs(np.cos(SIE[lat_label]*np.pi/180))
        SIE = SIE.sum(dim=(lat_label,lon_label))
        SIE = SIE.rename({sic_label: 'SIE'})
    else:
        SIE = grid_area.expand_dims(time=SIC['time']).where((SIC[sic_label] >= 0.15) & (SIC[sic_label] <= 1.0)).sum(dim=(lat_label,lon_label))
        SIE = SIE.rename({'cell_area': 'SIE'})
    SIE *= mfactor

    # average/select ensembles if specified
    if ensemble == 'average' or ensemble == 'ave' or ensemble == 'mean':
        SIE = SIE.mean(dim=ensemble_label)
    elif not ensemble == None:
        SIE = SIE.where(SIE[ensemble_label] == ensemble, drop=True)
        SIE = SIE.squeeze(dim=ensemble_label)

    # match time coordinate format with that of input and then return the dataset
    SIE['time'] = SIC['time']
    return SIE



def get_icearea (sic_dataset, grid_area_dataset=None, lat_bounds=None, lon_bounds=None, lat_label='lat', lon_label='lon', sic_label='SICN', flip_meridional=False, flip_zonal=False, ensemble=None, ensemble_label='ensemble', sic_factor=1, mfactor=1e-12):
    """
    Function that takes in a netcdf dataframe of gridded sea ice concentration (SIC) and outputs the corresponding sea ice area (in units of 10^6 km^2, unless specified by user w/ mfactor) 
    at each time step, conventionally defined as the integrated area of all grid cells having SIC > 0.15

    Args:
        sic_dataset:                an xarray or netcdf dataset containing a SIC variable as a function of coordinates latitude, longitude, and time.
        grid_area_dataset:          an xarray or netcdf dataset containing grid cell area as a function of coordinates latitude and longitude.
                                    Can be computed through the climate data operators command 'cdo gridarea'.
                                    If no grid is provided (i.e., left as None), then SIE will be calculated with a less accurate method specific to 1x1 grids.
        lat_bounds (tuple):         the lower and upper bound of latitudes in which SIC will be considered. Outside of the bounds, SIC is set to nan.
        lon_bounds (tuple):         similar to lat_bounds, except for longitude. Check whether your dataset inherently uses (0,360) or (-180,180).
        lat_label (string):         label of the latitude coordinate in the dataset (not the grid file). Often 'lat' or 'latitude'. 
        lon_label (string):         label of the longitude coordinate in the dataset (not the grid file). Often 'lon' or 'longitude'.
        sic_label (string):         label of the SIC data variable in the dataset. Often 'SICN' or 'siconc'.
        flip_meridional (bool):     if True, then sic INSIDE the LATITUDE bounds are set to nans, rather than outside. Defaults to False.
        flip_zonal (bool):          if True, then sic INSIDE the LONGITUDE bounds are set to nans, rather than outside. Defaults to False.
        ensemble (None, str, int):  selects the given ensemble and drops the rest. if 'ave' or 'mean', then averages over ensembles instead.
        ensemble_label (str):       label of the ensemble coordinate in the dataset. Often 'ensemble' or 'realization'.
        sic_factor (float):         multiplies all sic values by this factor. useful for converting from percents. default 1.
        mfactor (float):            multiplicative factor to multiple the final result by. Use 1e-12 if you want the outputted SIA to be in units of 10^6 sq km. If 1, then output is in units of sq m.

    Returns:
        (xarray dataset, optional)  the dataset of calculated SIA.
    """

    grid_area = grid_area_dataset
    SIC = sic_dataset
    if not type(grid_area_dataset) == type(None):
        if grid_area.dims == {'lon':360, 'lat':180}:
            grid_area = grid_area.rename({'lat':lat_label,'lon':lon_label})

    # apply latitude and longitude bounds
    SIC = _apply_bounds(SIC, lat_bounds, lon_bounds, lat_label, lon_label, flip_zonal, flip_meridional)

    # calculate SIA
    SIC[sic_label] *= sic_factor
    SIC = SIC.where((SIC[sic_label] <= 1.)&(SIC[sic_label] >= 0.),drop=True)
    if type(grid_area_dataset) == type(None):
        SIC[sic_label] = SIC[sic_label] * (111120**2)*np.abs(np.cos(SIC[lat_label]*np.pi/180))
    else:
        SIC[sic_label] = SIC[sic_label] * grid_area.expand_dims(time=SIC['time'])['cell_area']
    SIA = SIC.sum(dim=(lat_label,lon_label)).rename({sic_label: 'SIA'})
    SIA *= mfactor

    # average/select ensembles if specified
    if ensemble == 'average' or ensemble == 'ave' or ensemble == 'mean':
        SIA = SIA.mean(dim=ensemble_label)
    elif not ensemble == None:
        SIA = SIA.where(SIA[ensemble_label] == ensemble, drop=True)
        SIA = SIA.squeeze(dim=ensemble_label)

    # match time coordinate format with that of input and then return the dataset
    SIA['time'] = SIC['time']
    return SIA



def get_climatology (dataset, var, ref_period=None):
    """
    Function that takes in a netcdf dataset and calculates the MONTHLY climatology for a user-specified variable and reference period.

    Args:
        dataset (xarray dataset):   the dataset containing a time coordinate and the variable for which the climatology is to be calculated for.
        var (string):               the label of the variable for the climatology (e.g., "SIE")
        ref_period (tuple):         a tuple of two years defining the timespan (inclusive) over which the climatology will be calculated. If 'None', then entire time span is used.

    Returns:
        (xarray dataset):           The climatology dataset.
    """

    # include only the data within the reference period range
    if not ref_period == None:
        start_year, end_year = ref_period
        dataset = dataset.where((dataset['time.year'] >= start_year) & (dataset['time.year'] <= end_year))

    # calculate the climatology for each month
    climatology = dataset.groupby('time.month').mean(dim='time')
    climatology['overall_mean'] = dataset.mean(dim='time')[var]

    return climatology



def get_anomalies (dataset, var, ref_period=None, ref_dataset=None):
    """
    Function that takes in a netcdf dataset and return anomalies relative to a MONTHLY climatology (based on a mean) for a user-specified variable and reference period.

    Args:
        dataset (xarray):           the dataset containing a time coordinate and the variable from which anomalies will be calculated.
        var (string):               the label of the variable to calculate anomalies (e.g., "SIE")
        ref_period (tuple):         a tuple of two years defining the timespan (inclusive) over which the climatology will be calculated. If 'None', then entire time span is used.
        ref_dataset (xarray):       the dataset whose climatology will be used to calculate the anomalies. If None, then anomalies are
                                    calculated with respect to the input dataset.
        
    Returns:
        (dataset, same dims as input): dataset as anomalies, with respect to the climatology calculated over the reference period provided.
    """

    dataset = dataset.copy(deep=True)

    # climatology for all months
    if type(ref_dataset) == type(None):
        climatology = get_climatology(dataset=dataset, var=var, ref_period=ref_period)[var]
    else:
        climatology = get_climatology(dataset=ref_dataset, var=var, ref_period=ref_period)[var]

    # subtract climatology from dataset to produce anomalies
    dataset[var] = dataset[var].groupby('time.year').map(lambda x: x - climatology.where(climatology['month'] <= len(x), drop=True).values)
    return dataset



def get_trend (dataset):
    """
    Function that computes a least-squares linear fit for the given data variable, assuming the data are equidistant in time.
    Currently this function is only compatible with a 1-dimensional time series. (i.e., the only coordinate is time.)

    Args:
        dataset (xarray dataset): the dataset containing a time coordinate and the variable to be detrended.

    Returns:
        (1D list): the slope and the intercept of the fit.
    """

    N = len(dataset['time'])
    slope, intercept = np.polyfit(np.arange(0,N,1), dataset, deg=1)
    return [slope,intercept]



def remove_mean (dataset, var, ref_period=(0,9999)):
    """
    Function that removes the mean, producing anomalies, for each month separately. (I.e. the November trend is calculated only using November data points.)

    Args:
        dataset (xarray dataset):   the dataset containing a time coordinate and the variable to be detrended.
        var (string):               the label of the variable to calculate anomalies (e.g., "SIE")
        ref_period (tuple):         use if you want the trend to be calculated only using data from a range of years, e.g., between the years 2000 and 2020.
                                    Note that the trend is removed from the entire dataset, regardless of the reference period used for its calculation, UNLESS
                                    you use the 'Mitch' method, in which case only data in the reference period are detrended.

    Returns:
        (dataset, same dims as input): dataset with the mean removed.
    """

    dataset = dataset.copy(deep=True) # deep must be set to True

    for m in np.arange(1,13,1): # loop over each month

        # get subset of data for month 'm'
        month_indices = (dataset['time.month'] == m)
        subset_month = dataset.where(dataset['time.month'] == m,drop=True)
        subset_detrended = subset_month.copy(deep=True) # the dataset to be (but not yet) detrended then returned. 'deep' must be true.

        # another subset, this time narrowed down to the years over which the trend is to be computed.
        year_indices = (subset_month['time.year'] >= ref_period[0]) & (subset_month['time.year'] <= ref_period[1])
        subset_monthandyear = subset_month.where(year_indices, drop=True)

        # Remove the mean
        subset_detrended[var] = subset_month[var] - subset_monthandyear[var].mean(dim='time')
        dataset[var].loc[month_indices] = subset_detrended[var]

    return dataset



def remove_trend (dataset, var, ref_period=(0,9999)):
    """
    Function that removes linear trends, for each month separately. (I.e. the November trend is calculated only using November data points.)
    This function is only compatible with a 1-dimensional time series. (i.e., the only coordinate is time; no longitude or latitude)

    Args:
        dataset (xarray dataset):   the dataset containing a time coordinate and the variable to be detrended.
        var (string):               the label of the variable to calculate anomalies (e.g., "SIE")
        ref_period (tuple):         use if you want the trend to be calculated only using data from a range of years, e.g., between the years 2000 and 2020.
                                    Note that the trend is removed from the entire dataset, regardless of the reference period used for its calculation, UNLESS
                                    you use the 'Mitch' method, in which case only data in the reference period are detrended.

    Returns:
        (dataset, same dims as input): dataset with the trend removed.
    """

    dataset = dataset.copy(deep=True) # deep must be set to True

    for m in np.arange(1,13,1): # loop over each month

        # get subset of data for month 'm'
        month_indices = (dataset['time.month'] == m)
        subset_month = dataset.where(month_indices,drop=True)
        subset_detrended = subset_month.copy(deep=True) # the dataset to be (but not yet) detrended then returned. 'deep' must be true.

        # another subset, this time narrowed down to the years over which the trend is to be computed.
        year_indices = (subset_month['time.year'] >= ref_period[0]) & (subset_month['time.year'] <= ref_period[1])
        subset_monthandyear = subset_month.where(year_indices, drop=True)

        # Remove linear trend
        slope, intercept = get_trend(subset_monthandyear[var])
        subset_detrended[var] = subset_month[var] - (intercept + slope * np.arange(len(subset_month)))
        dataset[var].loc[month_indices] = subset_detrended[var]

    return dataset


