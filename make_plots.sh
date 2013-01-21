#!/bin/bash

set -x

R=f
GS=500
fill_value=-2e9

for gamma in 1.0 2.0 5.0 10.0 20.0 50.0
do

  for alpha in 0.0 1.0 2.0
  do
      for project in "jakobshavn" "79N" "helheim"
      do
      nc2cdo.py ${project}_alpha_${alpha}_gamma_${gamma}.nc
      nccopy ${project}_alpha_${alpha}_gamma_${gamma}.nc ${project}_alpha_${alpha}_gamma_${gamma}_masked.nc
      ncks -A -v thk ${project}_cresis_${GS}m.nc ${project}_alpha_${alpha}_gamma_${gamma}_masked.nc
      ncap2 -O -s "where(thk<=0.) {topg=$fill_value; divHU=$fill_value; divHU_cresis=$fill_value; divHU_umt=$fill_value; divHU_searise=$fill_value;}" ${project}_alpha_${alpha}_gamma_${gamma}_masked.nc ${project}_alpha_${alpha}_gamma_${gamma}_masked.nc
      ncatted -a _FillValue,topg,o,f,$fill_value -a _FillValue,divHU,o,f,$fill_value -a _FillValue,divHU_umt,o,f,$fill_value -a _FillValue,divHU_cresis,o,f,$fill_value -a _FillValue,divHU_searise,o,f,$fill_value ${project}_alpha_${alpha}_gamma_${gamma}_masked.nc
      ~/base/PyPISMTools/scripts/basemap-plot.py -p medium --geotiff_file MODIS${project}250m.tif --coastlines --map_resolution $R --colorbar_label -v topg --singlerow --colormap /Users/andy/base/PyPISMTools/colormaps/wiki-2.0.cpt -o ${project}_alpha_${alpha}_gamma_${gamma}_topg.png ${project}_alpha_${alpha}_gamma_${gamma}_masked.nc
      ~/base/PyPISMTools/scripts/basemap-plot.py -p medium --geotiff_file MODIS${project}250m.tif --coastlines --map_resolution $R --colorbar_label --bounds -100 100 -v divHU --singlerow --colormap RdBu_r -o ${project}_alpha_${alpha}_gamma_${gamma}_divHU.png ${project}_alpha_${alpha}_gamma_${gamma}_masked.nc
     done
  done
done

gamma=1.0
alpha=1.0

for project in "jakobshavn" "79N" "helheim"
do
    ~/base/PyPISMTools/scripts/basemap-plot.py -p medium --geotiff_file MODIS${project}250m.tif --coastlines --map_resolution $R --colorbar_label --bounds -100 100 -v divHU_cresis --singlerow --colormap RdBu_r -o ${project}_cresis_divHU.png ${project}_alpha_${alpha}_gamma_${gamma}_masked.nc
    ~/base/PyPISMTools/scripts/basemap-plot.py -p medium --geotiff_file MODIS${project}250m.tif --coastlines --map_resolution $R --colorbar_label --bounds -100 100 -v divHU_searise --singlerow --colormap RdBu_r -o ${project}_searise_divHU.png ${project}_alpha_${alpha}_gamma_${gamma}_masked.nc
    ~/base/PyPISMTools/scripts/basemap-plot.py -p medium --geotiff_file MODIS${project}250m.tif --coastlines --map_resolution $R --colorbar_label --bounds -100 100 -v divHU_umt --singlerow --colormap RdBu_r -o ${project}_umt_divHU.png ${project}_alpha_${alpha}_gamma_${gamma}_masked.nc
done

for project in "jakobshavn" "79N" "helheim"
do
    for dataset in "cresis" "umt" "searise_v1.1"
    do
        ~/base/PyPISMTools/scripts/basemap-plot.py -p medium --geotiff_file MODIS${project}250m.tif --coastlines --map_resolution $R --colorbar_label -v topg --singlerow --colormap /Users/andy/base/PyPISMTools/colormaps/wiki-2.0.cpt -o ${project}_${dataset}_${GS}m_topg.png ${project}_${dataset}_${GS}m.nc
    done
done
