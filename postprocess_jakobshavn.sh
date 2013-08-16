#!/bin/bash

set -x

R=l
GS=250
fill_value=-2e9

alpha=0.0

for gamma in 1.0 2.0 5.0 10.0 20.0 50.0 100.0 200.0 500.0 1000.0 2000.0 5000.0 10000.0 20000.0 50000.0
do
    for project in "jakobshavn"
    do
        for exp in  "nodhdt_nobmelt" "dhdt_nobmelt"
        do

	    dataset=flightlines
            basemap-plot.py -p medium --colorbar_label -v thk --singlerow --bounds -200 200 --colormap RdBu_r -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_${dataset}_thk_diff.pdf --obs_file ${project}_${dataset}_${GS}m.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
            extract-profile.py jakobshavn_xprofile_str_point.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc ${project}/profile_str_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
            extract-profile.py jakobshavn_xprofile_s8_point.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc ${project}/profile_s8_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
          done
     done
done

for project in "jakobshavn"
do
    for dataset in "cresis" "umt" "searise_v1.1"
    do
        # thk of input data sets
        basemap-plot.py -p medium --colorbar_label --bounds 0 2500 -v thk --singlerow --colormap /Users/andy/base/PyPISMTools/colormaps/sst.cpt -o ${project}/${project}_${dataset}_${GS}m_thk.pdf ${project}_${dataset}_${GS}m.nc
        # topg of input data sets
        basemap-plot.py -p medium --colorbar_label -v topg --singlerow --colormap /Users/andy/base/pypismtools/colormaps/wiki-2.0.cpt -o ${project}/${project}_${dataset}_${GS}m_topg.pdf ${project}_${dataset}_${GS}m.nc
        extract-profile.py jakobshavn_xprofile_str_point.shp ${project}_${dataset}_${GS}m.nc profile_str_${project}_${dataset}_${GS}m.nc
        extract-profile.py jakobshavn_xprofile_s8_point.shp ${project}_${dataset}_${GS}m.nc profile_s8_${project}_${dataset}_${GS}m.nc
    done
done
