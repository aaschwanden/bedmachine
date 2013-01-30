#!/bin/bash

set -x

R=l
GS=500
fill_value=-2e9


for gamma in 1.0 5.0 10.0 50.0 100.0
do
  for alpha in 5.0 10.0
  do
      for project in "jakobshavn"
      do
          for exp in "dhdt_bmelt" "dhdt_nobmelt" "nodhdt_bmelt" "nodhdt_nobmelt"
          do
              nc2cdo.py ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
              nccopy ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_masked.nc
              for var in "topg" "thk" "divHU"
              do
              gdal_translate -a_srs EPSG:3413 -of GTiff NETCDF:${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc:${var} ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_${var}.tif
              gdaldem hillshade -s 0.5 ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_$var.tif ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_${var}_hillshade.tif
              done
              ncks -A -v thk ${project}_cresis_${GS}m.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_masked.nc
              ncap2 -O -s "where(thk<=0.) {topg=$fill_value; divHU=$fill_value; divHU_cresis=$fill_value; divHU_umt=$fill_value; divHU_searise=$fill_value;}" ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_masked.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_masked.nc
              ncatted -a _FillValue,topg,o,f,$fill_value -a _FillValue,divHU,o,f,$fill_value -a _FillValue,divHU_umt,o,f,$fill_value -a _FillValue,divHU_cresis,o,f,$fill_value -a _FillValue,divHU_searise,o,f,$fill_value -a _FillValue,divU,o,f,$fill_value ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_masked.nc
      ## topg
              basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif --shape_file jakobshavn_terminus_2009 --map_resolution $R --colorbar_label -v topg --singlerow --colormap /Users/andy/base/PyPISMTools/colormaps/wiki-2.0.cpt -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_topg.png ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_masked.nc
      ## divHU
              basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif --map_resolution $R --colorbar_label -v divHU --singlerow --colormap RdBu_r  --bounds -100 100 -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_divHU.png ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_masked.nc
      ## thk
      # Don't plot 'thk' for masked files -- it's all the same because of the above masking operation
              basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif  --map_resolution $R --colorbar_label --bounds 0 2500 -v thk --singlerow --colormap Set1_r -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_thk.png ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
          done
     done
  done
done

gamma=1.0
alpha=0.0

for project in "jakobshavn"
do
    # thk flightlines
    basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif  --map_resolution $R --colorbar_label --bounds 0 2500 -v thk --singlerow --colormap Set1_r -o ${project}/${project}_${GS}m_flightlines_thk.png tmp_${project}_flightlines_${GS}m.nc
    # surf vels
    basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif  --map_resolution $R --colorbar_label -v magnitude --singlerow -o ${project}/${project}_${GS}m_speed.png ${project}_surf_vels_${GS}m.nc
    # divHU cresis
    basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif  --map_resolution $R --colorbar_label --bounds -100 100 -v divHU_cresis --singlerow --colormap RdBu_r -o ${project}/${project}_${GS}m_cresis_divHU.png ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_masked.nc
    # divHU searise
    basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif  --map_resolution $R --colorbar_label --bounds -100 100 -v divHU_searise --singlerow --colormap RdBu_r -o ${project}/${project}_${GS}m_searise_divHU.png ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_masked.nc
    # divHU UMT
    basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif  --map_resolution $R --colorbar_label --bounds -100 100 -v divHU_umt --singlerow --colormap RdBu_r -o ${project}/${project}_${GS}m_umt_divHU.png ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_masked.nc
done

for project in "jakobshavn"
do
    for dataset in "cresis" "umt" "searise_v1.1"
    do
        # thk of input data sets
        basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif  --map_resolution $R --colorbar_label --bounds 0 2500 -v thk --singlerow --colormap Set1_r -o ${project}/${project}_${dataset}_${GS}m_thk.png ${project}_${dataset}_${GS}m.nc
        # topg of input data sets
        basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif  --map_resolution $R --colorbar_label -v topg --singlerow --colormap /Users/andy/base/PyPISMTools/colormaps/wiki-2.0.cpt -o ${project}/${project}_${dataset}_${GS}m_topg.png ${project}_${dataset}_${GS}m.nc
    done
done

for gamma in 1.0 2.0 5.0 10.0 20.0 50.0
do
  for alpha in 0.0 1.0 2.0
  do
      for project in "jakobshavn"
      do
          for dataset in "cresis" "umt" "searise_v1.1"
          do
              # thk difference exp-obs
              basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif --map_resolution $R --colorbar_label -v thk --singlerow --bounds -200 200 --colormap RdBu_r -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_${dataset}_thk_diff.png --obs_file ${project}_${dataset}_${GS}m.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
              # topg difference exp-obs
              basemap-plot.py -p medium --geotiff_file ${project}_gimpdem_90m_hillshade.tif --map_resolution $R --colorbar_label -v topg --singlerow --bounds -200 200 --colormap RdBu_r -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_${dataset}_topg_diff.png --obs_file ${project}_${dataset}_${GS}m.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}.nc
           done
     done
  done
done

exit

