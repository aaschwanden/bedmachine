#!/bin/bash

set -x

for gamma in 1.0 2.0 5.0 10.0 20.0 50.0
do

  for alpha in 0.0 1.0 2.0
  do
      for project in "jakobshavn" "79N" "helheim"
      do
      nc2cdo.py ${project}_alpha_${alpha}_gamma_${gamma}.nc
      ~/base/PyPISMTools/scripts/basemap-plot.py -v topg --singlerow --colormap /Users/andy/base/PyPISMTools/colormaps/wiki-2.0.cpt -o ${project}_alpha_${alpha}_gamma_${gamma}_topg.png ${project}_alpha_${alpha}_gamma_${gamma}.nc
      ~/base/PyPISMTools/scripts/basemap-plot.py --bounds -100 100 -v divHU --singlerow --colormap RdBu_r -o ${project}_alpha_${alpha}_gamma_${gamma}_divHU.png ${project}_alpha_${alpha}_gamma_${gamma}.nc
     done
  done
done

gamma=1.0
alpha=1.0

for project in "jakobshavn" "79N" "helheim"
do
    ~/base/PyPISMTools/scripts/basemap-plot.py --bounds -100 100 -v divHU_cresis --singlerow --colormap RdBu_r -o ${project}_cresis_divHU.png ${project}_alpha_${alpha}_gamma_${gamma}.nc
    ~/base/PyPISMTools/scripts/basemap-plot.py --bounds -100 100 -v divHU_searise --singlerow --colormap RdBu_r -o ${project}_searise_divHU.png ${project}_alpha_${alpha}_gamma_${gamma}.nc
    ~/base/PyPISMTools/scripts/basemap-plot.py --bounds -100 100 -v divHU_umt --singlerow --colormap RdBu_r -o ${project}_umt_divHU.png ${project}_alpha_${alpha}_gamma_${gamma}.nc
done
