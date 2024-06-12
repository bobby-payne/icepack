#!/bin/bash

months=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")

# loop over init. month
for mo in ${months[@]}; do
	echo "Processing data for initialization month ${mo}"
	sleep 1
	# loop over realization
	for ens in $(seq 1 10); do

		# the data for the realization and initialization month
        	data=$(echo *${ens}i1p2f1_gn_*${mo}-*.nc)

        	cdo mergetime $data siconc_SImon_CanESM5_chfp3b-hindcast_s1980123024-r${ens}i1p2f1_gn_init1981${mo}.nc

	done
done
