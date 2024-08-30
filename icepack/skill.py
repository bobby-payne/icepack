### ---------------------------------------------
### Author: Robert Payne
### Contact: gpayne1654@uvic.ca
###
### This package makes extensive use of array indexing techniques, and the libraries numpy, xarray, and matplotlib. 
### It's recommended you familiarize yourself with these before working with this library.
### ---------------------------------------------

import numpy as np
import xarray as xr
from numba import njit



@njit
def _sample_withreplacement(dataset1, dataset2, seed=None):

    # apply seed, if specified
    if not seed == None:
        np.random.seed(seed)

    # check that datasets are of equal length
    N_data = len(dataset1)
    if not N_data == len(dataset2):
        raise IndexError("Datasets are not of equal length.")
    
    # generate N_data random indices, each between 0 and N_data, and select the corresponding subsets of data
    indices = np.empty(N_data, dtype=np.int32)
    for i in range(N_data):
        indices[i] = np.random.randint(N_data)
    # indices = np.random.randint(N_data, size=N_data)
    datasubset1 = dataset1[indices]
    datasubset2 = dataset2[indices]

    return [datasubset1, datasubset2]



def get_skillmatrix (data_obs, data_initialized, daterange=None, significance_test=False, significance_level=95, Nb=1000, seed=None):

    # empty skill matrix, to be filled
    skill_matrix = np.full((12,12),np.nan)
    if significance_test:
        significance_matrix = np.full((12,12),False)
    if not type(daterange) == type(None):
        y0,y1 = daterange
    if 'lead' in data_initialized[0].data_vars:
        data_initialized = [data.drop_vars('lead') for data in data_initialized]

    # loop over initialized data
    for i,data_init in enumerate(data_initialized):
        # loop over target month
        for tm in np.arange(1,13):

            # data for only that target month
            data_init_subset = data_init.where(data_init['time.month']==tm, drop=True)
            data_obs_subset = data_obs.where(data_obs['time.month']==tm, drop=True)

            # apply daterange
            if not type(daterange) == type(None):
                data_init_subset = data_init_subset.where((data_init_subset['time.year']>=y0)&(data_init_subset['time.year']<=y1),drop=True)
                data_obs_subset = data_obs_subset.where((data_obs_subset['time.year']>=y0)&(data_obs_subset['time.year']<=y1),drop=True)
                data_init_subset['time'] = data_obs_subset['time']
            else:
            # align along the time axis
                data_init_subset, data_obs_subset = xr.align(data_init_subset, data_obs_subset)
        
            # compute skill
            skill = xr.corr(data_obs_subset['SIE'], data_init_subset['SIE'], dim='time').values.item()
            skill_matrix[i][tm-1] = skill

            # compute significance of the skill, if significance_test is True
            if significance_test:

                # sample with replacement Nb times
                data_init_subset_array = data_init_subset['SIE'].values
                data_obs_subset_array = data_obs_subset['SIE'].values
                bootstrapped_samples = [_sample_withreplacement(data_init_subset_array, data_obs_subset_array, seed=seed) for _ in range(Nb)]
                bootstrapped_correlations = np.array([np.corrcoef(sample[0], sample[1])[0][1] for sample in bootstrapped_samples])

                # calculate the proportion of correlations greater than (less than) zero, and use to determine if significant.
                proportion = np.sum(np.sign(bootstrapped_correlations) == np.sign(skill))/Nb
                if proportion >= significance_level/100.:
                    significance_matrix[i][tm-1] = True

    if significance_test:
        return [skill_matrix, significance_matrix]
    else:
        return skill_matrix


