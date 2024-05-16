### ---------------------------------------------
### Author: Robert Payne
### Contact: gpayne1654@uvic.ca
###
### This package makes extensive use of array indexing techniques, and the libraries numpy, pandas, and xarray. 
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

