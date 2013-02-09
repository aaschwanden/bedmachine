#!/bin/bash

set -x -e

for VYEARS in "2008_2009" "2007_2008" "2006_2007"; do
# remap surface velocities, select area frist to speed
# things up a bit
    VELIN_FILE=surf_vels_500m_${VYEARS}.nc
    VELOUT_FILE=${PROJECT}_surf_vels_${VYEARS}_${GS}m.nc
    #cdo setmisstoc,0 $VELIN_FILE tmp_${VELIN_FILE}
    #ncap2 -O -s "where(ue>100) {us=-2e9; vs=-2e9; magnitude=-2e9;};" tmp_${VELIN_FILE} tmp_${VELIN_FILE}
    #fill_missing.py -v magnitude,us,vs -e 1 -f tmp_${VELIN_FILE} -o tmp2_${VELIN_FILE}
    if [ [$NN == 1] ] ; then
        cdo remapbil,$FL_FILE_NC -sellonlatbox,$LON_MIN,$LON_MAX,$LAT_MIN,$LAT_MAX tmp_2$VELIN_FILE $VELOUT_FILE
    else
        cdo -P $NN remapbil,$FL_FILE_NC  -sellonlatbox,$LON_MIN,$LON_MAX,$LAT_MIN,$LAT_MAX tmp2_$VELIN_FILE $VELOUT_FILE
    fi
    ncks -A -v x,y,mapping $FL_FILE_NC $VELOUT_FILE
    ncatted -a grid_mapping,us,o,c,"mapping" -a grid_mapping,vs,o,c,"mapping" -a grid_mapping,magnitude,o,c,"mapping" $VELOUT_FILE
done

for VYEARS in "2008_2009" "2007_2008" "2006_2007"; do
# remap surface velocities, select area frist to speed
# things up a bit
    VELIN_FILE=surf_vels_500m_${VYEARS}.nc
    VELOUT_FILE=${PROJECT}_surf_vels_${VYEARS}_${GS}m.nc
    python scripts/interpolate-velocities.py -o tmp3_$VELIN_FILE $VELIN_FILE
    #cdo setmisstoc,0 $VELIN_FILE tmp_${VELIN_FILE}
    #ncap2 -O -s "where(ue>100) {us=-2e9; vs=-2e9; magnitude=-2e9;};" tmp_${VELIN_FILE} tmp_${VELIN_FILE}
    #fill_missing.py -v magnitude,us,vs -e 1 -f tmp_${VELIN_FILE} -o tmp2_${VELIN_FILE}
    if [ [$NN == 1] ] ; then
        cdo remapbil,$FL_FILE_NC -sellonlatbox,$LON_MIN,$LON_MAX,$LAT_MIN,$LAT_MAX tmp3_$VELIN_FILE int_$VELOUT_FILE
    else
        cdo -P $NN remapbil,$FL_FILE_NC  -sellonlatbox,$LON_MIN,$LON_MAX,$LAT_MIN,$LAT_MAX tmp3_$VELIN_FILE int_$VELOUT_FILE
    fi
    ncks -A -v x,y,mapping $FL_FILE_NC int_$VELOUT_FILE
    ncatted -a grid_mapping,us,o,c,"mapping" -a grid_mapping,vs,o,c,"mapping" -a grid_mapping,magnitude,o,c,"mapping" int_$VELOUT_FILE
done
