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
import pandas as pd
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
import netCDF4, h5netcdf, dask
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib import cm, colors, colormaps
from scipy.io import loadmat

def _add_gridlines (ax, color='black', alpha=1.):
    """
    For internal use only. Adds lines to `ax` to denote the boundaries of the regions.
    As far as I'm concerned, this function is black magic.
    """
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,
                  linewidth=1, color=color, alpha=alpha, linestyle='--',dms=True)
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.xlines = True
    gl.ylines = False
    gl.xlocator = mticker.FixedLocator([20,90,160,-130,-60])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 15, 'color': 'grey'}
    gl.xlabel_style = {'color': 'red', 'weight': 'bold'}



def plot_fields (var_data, var_names, cbarlabel="Sea Ice Concentration", clevels=None, clevelcolor='grey', figsize=None, region_lines=False, return_axes=False):
    """
    Plots two fields around the Antarctic, one as a colormesh and another as contours (e.g., for sic and psl).

    Args:
        var_data (list):        var_data must be a list of TWO 2D arrays. Each 2D array contains the data sets to be plotted.
                                The first dataset will become the colormesh, and the second
                                will become the contours. Each dataset MUST consist of a single gridded variable as a function of longitude and latitude coordinates,
                                but NO time coordinate (i.e., choose a specific time before using this function).
        var_names (list):       the names of the variables to be plotted, as in the xarray datasets (e.g., ['sic', 'psl']).
        cbarlabel (str):        Label for the colormesh colour bar
        clevels (1D array):     A 1D array defining the contour levels for the contour plot. If None, then will be assigned automatically.
        clevelcolor (str):      the colour of the contour lines in the plot. Defaults to grey.
        figsize (tuple):        The size of the figure, as in matplotlib (e.g., (8,6)). If None, then will be assigned automatically.
        region_lines (bool):    Whether or not you want dashed lines added to denote the borders of the different Antarctic regions. Defaults to false.
        return_axes (bool):     Whether or not you want this function to return the matplotlib axes list for the plot. Defaults to False.
    """

    # process input data
    var1_matrix, var2_matrix = var_data
    var1_name, var2_name = var_names
    nrows, ncols = [len(var1_matrix), len(var1_matrix[0])]
    if not ((len(var1_matrix) == len(var2_matrix)) and (len(var1_matrix[0]) == len(var2_matrix[0]))):
        raise ValueError("Input matrices are not of same shape.")
        
    # circular plots rather than square
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], .5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    # plot parameters/design
    cmap = colors.ListedColormap(loadmat('./cmaps/cmap_jet3.mat')['cmap'], name='jet3')
    cmap.set_bad(color='grey', alpha=1.0)  # Specify the color for NaN values
    cticks = np.arange(-1,1.2,.2)
    projection = ccrs.Orthographic(central_latitude=-90., central_longitude=0.0)
    transform = ccrs.PlateCarree()                                              
    domain = [-180,180,-90,-47.5]

    # plot
    fig, ax = plt.subplots(nrows=nrows,ncols=ncols,figsize=figsize,facecolor='white',dpi=200,subplot_kw={'projection':projection})
    for i in range(nrows):
        for j in range(ncols):

            # get data from matrix
            var1 = var1_matrix[i][j]
            var2 = var2_matrix[i][j]

            # determine axis to plot on
            if nrows == 1 and ncols == 1:
                ax_ij = ax
            elif nrows == 1:
                ax_ij = ax[j]
            elif ncols == 1:
                ax_ij = ax[i]
            else:
                ax_ij = ax[i][j]

            # plot data
            ax_ij.set_extent(domain,transform)
            ax_ij.coastlines()
            if region_lines == True:
                _add_gridlines(ax_ij)
            var1_plot = var1[var1_name].plot(ax=ax_ij,transform=transform,cmap=cmap,vmin=-1,vmax=1,cbar_kwargs={"label":cbarlabel,"ticks":cticks})
            var2_plot = var2[var2_name].plot.contour(ax=ax_ij,transform=transform,colors=clevelcolor,levels=clevels,linewidths=.5)
            plt.clabel(var2_plot, inline=True, fontsize=3)

    plt.show()
    plt.tight_layout()

    if return_axes:
        return ax
