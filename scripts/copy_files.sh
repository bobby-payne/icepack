#!/bin/bash

ens=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10')
dest=/space/hall5/sitestore/eccc/crd/ccrn/users/rrp000/data/model/CanESM5/nudge+init/winds_201601-201612/psl_g/

for e in "${ens[@]}"
do
	if [ ${#e} -eq 1 ]; then
		cp /space/hall6/sitestore/eccc/crd/ccrn/users/rms101/c3bnw2-201601e0${e}/data/nc_output/CMIP6/CCCma/CCCma/CanESM5-c3bnw2-201601e0${e}/dcppA-hindcast/s2015-r${e}i1p2f1/Amon/psl/gn/v20190429/psl_Amon_CanESM5-c3bnw2-201601e0${e}_dcppA-hindcast_s2015-r${e}i1p2f1_gn_201601-201612.nc ${dest}
	else
		cp /space/hall6/sitestore/eccc/crd/ccrn/users/rms101/c3bnw2-201601e${e}/data/nc_output/CMIP6/CCCma/CCCma/CanESM5-c3bnw2-201601e${e}/dcppA-hindcast/s2015-r${e}i1p2f1/Amon/psl/gn/v20190429/psl_Amon_CanESM5-c3bnw2-201601e${e}_dcppA-hindcast_s2015-r${e}i1p2f1_gn_201601-201612.nc ${dest}
	fi
done

