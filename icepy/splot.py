### ---------------------------------------------
### Author: Robert Payne
### Contact: gpayne1654@uvic.ca
###
### This package makes extensive use of array indexing techniques, and the libraries numpy, pandas, xarray, and matplotlib. 
### It's recommended you familiarize yourself with these before working with this library.
### ---------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import xarray as xr
import cftime
import pandas as pd
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
import netCDF4, h5netcdf, dask
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib import cm, colors, colormaps
from scipy.io import loadmat



def set_extent(ax, domain = [-180,180,-90,-47.5], transform=ccrs.PlateCarree()):
    """
    Sets the latitudinal and longitudinal extent of the plot.
    """
    ax.set_extent(domain,transform)



def add_sic(ax, data, anom=True, month=None, year=None, transform=ccrs.PlateCarree()):
    """
    Adds a sea ice concentration heatmap to the given axis.
    """

    # select given year and month (otherwise, assumes user already pre-selected)
    if not year == None:
        data = data.where(data['time.year']==year,drop=True)
    if not month == None:
        try:
            data = data.where(data['time.month']==month,drop=True)
        except:
            data = data.where(data['month']==month,drop=True)

    # Get sic data from input dataset
    # Note to self: figure out what this code block does and why it's needed
    if 'siconc' in list(data.variables):
        sicname = 'siconc'
    else:
        sicname = 'SICN'
    try:
        data[sicname] = data[sicname][list(data[sicname].dims).index('time')] 
    except:
        data[sicname] = data[sicname][list(data[sicname].dims).index('month')]
    
    # colour map
    if anom:
        cmap = colors.ListedColormap(loadmat('./cmaps/cmap_jet3.mat')['cmap'], name='jet3')
        vmin,vmax = [-1,1]
        cbarlabel = "SIC Anomaly"
    else:
        cmap = colors.ListedColormap(loadmat('./cmaps/cmap_jet3_pos.mat')['cmap'], name='jet3')
        vmin,vmax = [0,1]
        cbarlabel = "SIC"
    cmap.set_bad(color='lightgrey', alpha=1)  # Specify the color for NaN values

    # plot
    sic_plot = data[sicname].plot(ax=ax,transform=transform,cmap=cmap,add_colorbar=False)
    return sic_plot


def add_ice_edge(ax, data, month=None, year=None, clevel=0.15, ccol='black', cls='--', clw=0.5, transform=ccrs.PlateCarree()):
    """
    Adds a sea ice edge contour to the given axis.
    """

    # select given year and month (otherwise, assumes user already pre-selected)
    if not year == None:
        data = data.where(data['time.year']==year,drop=True)
    if not month == None:
        try:
            data = data.where(data['time.month']==month,drop=True)
        except:
            data = data.where(data['month']==month,drop=True)

    # Get sic data from input dataset
    # Note to self: figure out what this code block does and why it's needed
    if 'siconc' in list(data.variables):
        sicname = 'siconc'
    else:
        sicname = 'SICN'
    try:
        data[sicname] = data[sicname][list(data[sicname].dims).index('time')] 
    except:
        data[sicname] = data[sicname][list(data[sicname].dims).index('month')]
    
    # plot
    ice_edge_plot = data[sicname].plot.contour(ax=ax,transform=transform,colors='black',levels=[clevel],linewidths=clw,linestyles=cls)
    return ice_edge_plot


def add_psl(ax, data, month=None, year=None, clevels=np.arange(950, 1200, 5), ccol='cyan', clw=0.5, transform=ccrs.PlateCarree()):
    """
    Adds sea level pressure contours to the given axis.
    """
    
    # select given year and month (otherwise, assumes user already pre-selected)
    if not year == None:
        data = data.where(data['time.year']==year,drop=True)
    if not month == None:
        try:
            data = data.where(data['time.month']==month,drop=True)
        except:
            data = data.where(data['month']==month,drop=True)
    try:
        data['psl'] = data['psl'][list(data.dims).index('time')]
    except:
        data['psl'] = data['psl'][list(data.dims).index('month')]

    # plot
    psl_plot = data['psl'].plot.contour(ax=ax,transform=transform,colors=ccol,levels= clevels,linewidths=clw)
    plt.clabel(psl_plot, inline=True, fontsize=4, colors='black')
    return psl_plot
