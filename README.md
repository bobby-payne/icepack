IcePy is a small package I developed for my research on sea ice at the University of Victoria, BC, Canada. Its functionality is makes heavy use of `numpy` and `xarray`.

Note that the folder `old_code` is for code I wrote for my undergraduate thesis. It is bulky, slow, and not very well documented -- I heavily advise against its use. The directory will be removed in the future, though for now it serves as an occassional reference and refresher on what I had done previously.

# How to convert SIC to SIE

To start, you should have a netcdf file with the sea ice concentration (SIC) data as a function of longitude and latitude, at the very least. For the sake of this README, I'll assume the file is called `sic_data.nc` and is located at `~\project\data\sic\sic_data.nc`. 

You will also need a netcdf file containing the grid cell areas. This is most easily obtained using the `cdo gridarea` command in shell. I'll call it `gridarea.nc` and store it in `~\project\data\sic\gridarea.nc`.

There are two commands you will want to import from my icepy package: `format_time_coord` and `sic_to_sie`. You will also need xarray. If your working directory has the icepy package in it, you can import these via
```
import xarray
from icepy.analysis import format_time_coord, sic_to_sie
```

The following code then computes SIC
```
# the locations of the data and grid area file
data_path = "~\project\data\sic\sic_data.nc"
grid_path = "~\project\data\sic\gridarea.nc"

# open data and convert the time coordinate into a cftime type.
# in this case, the time series starts in January 1980 and ends in December 2020.
sic = format_time_coord(xr.open_dataset(data_path),'1980-01','2020-12',freq='M')
grid = xr.open_dataset(grid_path)

# calculate SIE from SIC
# the lat_bounds in this case specify to only use the southern hemisephere.
# sic_label, lat_label, lon_label tell the function what the variables are named in the netcdf file.
sie = sic_to_sie(sic, grid, lat_bounds=(-90,0), sic_label='SICN', lat_label='lat', lon_label='lon')
```

Note that `sic_to_sie` (and other functions in this package) typically have many arguments that can be passed through them such that it is as flexible as possible (e.g., if you want the SIE to be averaged over multiple ensembles). Please refer to the docstrings of each function for a detailed description of the arguments that may be passed. 
