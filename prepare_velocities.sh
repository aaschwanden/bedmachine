#!/bin/bash

set -x -e

for VYEARS in "2008_2009" "2007_2008" "2006_2007"; do
# remap surface velocities, select area frist to speed
# things up a bit
    VELIN_FILE=surf_vels_500m_${VYEARS}.nc
    VELOUT_FILE=${PROJECT}_surf_vels_${VYEARS}_${GS}m.nc
    ncks -O -d x,$X_MIN,$X_MAX -d y,$Y_MIN,$Y_MAX $VELIN_FILE tmp_${VELIN_FILE}
    cdo -O setmisstoc,0 tmp_${VELIN_FILE} tmp2_${VELIN_FILE}
    ncap2 -O -s "where(ue>20) {us=-2e9; vs=-2e9; magnitude=-2e9;};" tmp2_${VELIN_FILE} tmp2_${VELIN_FILE}
    fill_missing.py -v magnitude,us,vs -e 1e-1 -f tmp2_${VELIN_FILE} -o tmp3_${VELIN_FILE}
    if [ [$NN == 1] ] ; then
        REMAP_EXTRAPOLATE=on cdo remapbil,$FL_FILE_NC tmp3_$VELIN_FILE $VELOUT_FILE
    else
        REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,$FL_FILE_NC tmp3_$VELIN_FILE $VELOUT_FILE
    fi
    ncks -A -v x,y,mapping $FL_FILE_NC $VELOUT_FILE
    ncatted -a grid_mapping,us,o,c,"mapping" -a grid_mapping,vs,o,c,"mapping" -a grid_mapping,magnitude,o,c,"mapping" $VELOUT_FILE
done

