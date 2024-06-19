### ---------------------------------------------
### Author: Robert Payne
### Contact: gpayne1654@uvic.ca
###
### This package makes extensive use of array indexing techniques, and the libraries numpy, xarray, and matplotlib. 
### It's recommended you familiarize yourself with these before working with this library.
### ---------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib import colors, colormaps
from scipy.io import loadmat

# Custom package
from icepy.analysis import select_date



def set_extent(ax, domain = [-180,180,-90,-47.5], transform=ccrs.PlateCarree()):
    """
    Sets the latitudinal and longitudinal extent of the plot.
    """
    ax.set_extent(domain,transform)



def add_sic(ax, data, anom=True, year=None, month=None, day=None, cmap=None, transform=ccrs.PlateCarree()):
    """
    Adds a sea ice concentration heatmap to the given axis.
    """

    # select given year and month (otherwise, assumes user already pre-selected)
    data = select_date(data,year=year,month=month,day=day,drop=True)

    # Get sic data from input dataset
    # Note to self: figure out what this code block does and why it's needed
    if 'siconc' in list(data.variables):
        sicname = 'siconc'
    elif 'SICN' in list(data.variables):
        sicname = 'SICN'
    else:
        sicname = 'sicn'
    
    # colour map
    if cmap == None:
        if anom:
            cmap = colormaps['seismic'](np.linspace(0,1,79))
            cmap = np.insert(cmap, 39, np.array([1.,1.,1.,1.]),axis=0) # add another white in the middle
            cmap = colors.ListedColormap(cmap)
            # cmap = colors.ListedColormap(loadmat('./cmaps/cmap_jet3.mat')['cmap'], name='jet3')
        else:
            cmap = colormaps['OrRd'](np.linspace(0,1,40))
            cmap[0] = np.array([1.,1.,1.,1.]) # zero is always white
            cmap = colors.ListedColormap(cmap)
            # cmap = colors.ListedColormap(loadmat('./cmaps/cmap_jet3_pos.mat')['cmap'], name='jet3')
    cmap.set_bad(color='lightgrey', alpha=1)  # Specify the color for NaN values
    cmap.set_extremes(over='white')

    # plot
    sic_plot = data[sicname].plot(ax=ax,transform=transform,cmap=cmap,add_colorbar=False)
    return sic_plot


def add_ice_edge(ax, data, year=None, month=None, day=None, clevel=0.15, ccol='black', cls='--', clw=0.5, transform=ccrs.PlateCarree()):
    """
    Adds a sea ice edge contour to the given axis.
    """

    # select given year and month (otherwise, assumes user already pre-selected)
    data = select_date(data,year=year,month=month,day=day,drop=True)

    # Get sic data from input dataset
    # Note to self: figure out what this code block does and why it's needed
    if 'siconc' in list(data.variables):
        sicname = 'siconc'
    elif 'SICN' in list(data.variables):
        sicname = 'SICN'
    else:
        sicname = 'sicn'
    try:
        data[sicname] = data[sicname][list(data[sicname].dims).index('time')] 
    except:
        data[sicname] = data[sicname][list(data[sicname].dims).index('month')]
    
    # plot
    ice_edge_plot = data[sicname].plot.contour(ax=ax,transform=transform,colors=ccol,levels=[clevel],linewidths=clw,linestyles=cls)
    return ice_edge_plot



def add_psl(ax, data, year=None, month=None, day=None, clevels=np.arange(950, 1200, 5), ccol='cyan', clw=0.5, transform=ccrs.PlateCarree()):
    """
    Adds sea level pressure contours to the given axis.
    """
    
    # select given year and month (otherwise, assumes user already pre-selected)
    data = select_date(data,year=year,month=month,day=day,drop=True)
    try:
        data['psl'] = data['psl'][list(data.dims).index('time')]
    except:
        data['psl'] = data['psl'][list(data.dims).index('month')]

    # plot
    try:
        data = data.assign_coords(lon=(((data.lon + 180) % 360) - 180)).sortby('lon')
    except:
        data = data.assign_coords(longitude=(((data.longitude + 180) % 360) - 180)).sortby('longitude')
    psl_plot = data['psl'].plot.contour(ax=ax,transform=transform,colors=ccol,levels= clevels,linewidths=clw)
    plt.clabel(psl_plot, inline=True, fontsize=4, colors='black', use_clabeltext=True)
    return psl_plot
