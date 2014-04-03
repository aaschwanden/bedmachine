#!/bin/bash

set -x -e
YEAR=1985
VELIN_FILE=surf_vels_${YEAR}_utm.nc
VELOUT_FILE=${PROJECT}_surf_vels_${YEAR}_${GS}m.nc
if [[ $NN == 1 ]] ; then
    REMAP_EXTRAPOLATE=on cdo remapbil,$FL_FILE_NC $VELIN_FILE $VELOUT_FILE
else
    REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,$FL_FILE_NC $VELIN_FILE $VELOUT_FILE
fi
ncks -A -v x,y,mapping $FL_FILE_NC $VELOUT_FILE
ncatted -a grid_mapping,us,o,c,"mapping" -a grid_mapping,vs,o,c,"mapping" -a grid_mapping,magnitude,o,c,"mapping" $VELOUT_FILE

