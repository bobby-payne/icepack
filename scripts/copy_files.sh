#!/bin/bash

ens=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10')
dest=~/sea_ice/data/model/CanESM5/nudge+init/day/winds+temp_20160101-20161231/psl/g/

for e in "${ens[@]}"
do
	if [ ${#e} -eq 1 ]; then
		cp /space/hall6/sitestore/eccc/crd/ccrn/users/rms101/c3bnwt2-201601e0${e}/data/nc_output/CMIP6/CCCma/CCCma/CanESM5-c3bnwt2-201601e0${e}/dcppA-hindcast/s2015-r${e}i1p2f1/day/psl/gn/v20190429/psl_day_CanESM5-c3bnwt2-201601e0${e}_dcppA-hindcast_s2015-r${e}i1p2f1_gn_20160101-20161231.nc ${dest}
	else
		cp /space/hall6/sitestore/eccc/crd/ccrn/users/rms101/c3bnwt2-201601e${e}/data/nc_output/CMIP6/CCCma/CCCma/CanESM5-c3bnwt2-201601e${e}/dcppA-hindcast/s2015-r${e}i1p2f1/day/psl/gn/v20190429/psl_day_CanESM5-c3bnwt2-201601e${e}_dcppA-hindcast_s2015-r${e}i1p2f1_gn_20160101-20161231.nc ${dest}
	fi
done

