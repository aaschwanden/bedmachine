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
        for exp in  "nodhdt_nobmelt"
        do
            nc2cdo.py ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
            # bed elevation
            basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif --map_resolution $R --colorbar_label -v topg --singlerow --colormap ~/base/PyPISMTools/colormaps/wiki-2.0.cpt -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_topg.png ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
            # thickness  
            basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif  --map_resolution $R --colorbar_label --bounds 0 2500 -v thk --singlerow --colormap ~/base/PyPISMTools/colormaps/sst.cpt -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_thk.png ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
            # thk difference
            basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif --map_resolution $R --colorbar_label -v thk --singlerow --colormap RdBu_r  --bounds -100 100 --obs_file ${project}_flightlines_${GS}m.nc -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_thkdiff.png ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
            # divHU
            basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif --map_resolution $R --colorbar_label -v divHU --singlerow --colormap RdBu_r  --bounds -100 100 -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_divHU.png ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
            extract-flowline.py jakobshavn_xprofile_str_point.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc ${project}/profile_str_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
            extract-flowline.py jakobshavn_xprofile_s8_point.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc ${project}/profile_s8_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
          done
     done
done

for project in "jakobshavn"
do
    for dataset in "cresis" "umt" "searise_v1.1"
    do
        # thk of input data sets
        basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif  --map_resolution $R --colorbar_label --bounds 0 2500 -v thk --singlerow --colormap /Users/andy/base/PyPISMTools/colormaps/sst.cpt -o ${project}/${project}_${dataset}_${GS}m_thk.png ${project}_${dataset}_${GS}m.nc
        # topg of input data sets
        basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif  --map_resolution $R --colorbar_label -v topg --singlerow --colormap /Users/andy/base/PyPISMTools/colormaps/wiki-2.0.cpt -o ${project}/${project}_${dataset}_${GS}m_topg.png ${project}_${dataset}_${GS}m.nc
        extract-flowline.py jakobshavn_xprofile_str_point.shp ${project}_${dataset}_${GS}m.nc profile_str_${project}_${dataset}_${GS}m.nc
        extract-flowline.py jakobshavn_xprofile_s8_point.shp ${project}_${dataset}_${GS}m.nc profile_s8_${project}_${dataset}_${GS}m.nc
    done
done
