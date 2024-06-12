import xarray as xr
import numpy as np
monthstr = ['01','02','03','04','05','06','07','08','09','10','11','12']
for m in monthstr:
	print(f"Working on init month {m}/12")
	data_all = [xr.open_dataset(f"~/sea_ice/data/model/CanESM5/init/sic/1x1/siconc_SImon_CanESM5_chfp3b-hindcast_s1980123024-r{e}i1p2f1_gn_init1981{m}.nc") for e in np.arange(1,11,1)]
	data = xr.concat(data_all,dim='ensemble')
	data.to_netcdf(f"~/sea_ice/data/model/CanESM5/init/sic/1x1/siconc_SImon_CanESM5_chfp3b-hindcast_s1980123024-i1p2f1_gn_init1981{m}.nc")
print("Done.")
