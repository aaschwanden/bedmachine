#!/bin/bash

# Copyright (C) 2013 Andy Aschwanden
#
# Description:
# script will download and process data sets, including
# querying ice thickness data base, gimp dem, cresis,
# SeaRISE and UMT data sets.
# surf_vels.nc containing surface velocities from SAR
# provided by Ian Joughin is also needed.
#
# Usage:
# This script is not stand-alone, but expected to be sourced
# from another bash script containing area extends and other
# parameters needed.
#
# You need to open a port in a different shell first:
# ssh -L5433:icebridge.sr.unh.edu:5432 bmachuser@icebridge.sr.unh.edu
# and public key exchange is required prior to usage.
#
# Returns:
# - all input files needed to run mass conserving bed estimator
# code provided by Jesse Johnson.

# FL_FILE_TXT=${PROJECT}_flightlines.txt
FL_FILE_TXT=${PROJECT}_cresis_flightlines_${YEARA}-${YEARE}.csv

# Well this takes a while...
# But we leave it at that for now. In the longer run, we need a python 
# script that downloads the data and puts it into a more appriate format
# than ASCII.

# python scripts/general_query.py -table cresis_gr -fields "thick,quality,frame" -year_range $YEARA $YEARE \
#    -epsg $EPSG -and_clause "(quality>-1 or quality<5) and (thick>-9999)" -box $LON_MIN $LON_MAX $LAT_MIN $LAT_MAX \
#    -mod_val $MOD_VAL -mod_field $MOD_FIELD > $FL_FILE_TXT

FL_FILE_NC=${PROJECT}_flightlines_${GS}m.nc
python scripts/resample-cresis-data.py -g $GS --bounds $X_MIN $X_MAX $Y_MIN $Y_MAX \
    -n $NN $FL_FILE_TXT tmp_$FL_FILE_NC
nc2cdo.py --srs '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m' tmp_$FL_FILE_NC

# nc2cdo.py is from pism/util/
# it adds lat/lon, but also the 4 grid corners of each cell, needed for
# conservative remapping via CDO.

# TODO figure out if we really need this.
cdo setmisstoc,-9999. -selvar,thk tmp_$FL_FILE_NC $FL_FILE_NC
ncatted -a _FillValue,,d,, $FL_FILE_NC
ncks -A -v thk -x tmp_$FL_FILE_NC $FL_FILE_NC
nc2cdo.py $FL_FILE_NC


WARPOPTIONS="-overwrite -multi -r bilinear -te $X_MIN $Y_MIN $X_MAX $Y_MAX -tr $GS $GS -t_srs EPSG:$EPSG"

SPOT_FILE_IN=jakobshavn_spot_dem_diff_clean.nc
SPOT_FILE_NC=${PROJECT}_dhdt_${GS}m.nc
if [ [$NN == 1] ] ; then
  REMAP_EXTRAPOLATE=on cdo remapbil,$FL_FILE_NC $SPOT_FILE_IN $SPOT_FILE_NC
else
  REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,$FL_FILE_NC $SPOT_FILE_IN $SPOT_FILE_NC
fi
ncrename -v Band1,dhdt $SPOT_FILE_NC
ncks -A -v x,y,mapping $FL_FILE_NC $SPOT_FILE_NC

BMELT_FILE_IN=g1km_0_CLRUN_bmelt.nc
BMELT_FILE_NC=${PROJECT}_bmelt_${GS}m.nc
if [ [$NN == 1] ] ; then
  REMAP_EXTRAPOLATE=on cdo remapbil,$FL_FILE_NC $BMELT_FILE_IN $BMELT_FILE_NC
else
  REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,$FL_FILE_NC $BMELT_FILE_IN $BMELT_FILE_NC
fi
ncks -A -v x,y,mapping $FL_FILE_NC $BMELT_FILE_NC


# CReSIS data set
CRESIS=${PROJECTC}_${CRESIS_YEARS}
CRESIS_FILE_ZIP=${CRESIS}_Composite.zip
CRESIS_FILE_NC=${PROJECT}_cresis_${GS}m.nc
wget -nc --no-check-certificate https://data.cresis.ku.edu/data/grids/$CRESIS_FILE_ZIP
unzip -o $CRESIS_FILE_ZIP
gdalwarp $WARPOPTIONS -of netCDF ${CRESIS}_Composite/grids/${PROJECT}_${CRESIS_YEARS}_composite_thickness.txt thk_$CRESIS_FILE_NC
gdalwarp $WARPOPTIONS -of netCDF ${CRESIS}_Composite/grids/${PROJECT}_${CRESIS_YEARS}_composite_surface.txt usurf_$CRESIS_FILE_NC
ncrename -O -v Band1,usurf usurf_$CRESIS_FILE_NC
gdalwarp $WARPOPTIONS -of netCDF ${CRESIS}_Composite/grids/${PROJECT}_${CRESIS_YEARS}_composite_bottom.txt topg_$CRESIS_FILE_NC
ncrename -O -v Band1,topg topg_$CRESIS_FILE_NC
nccopy thk_$CRESIS_FILE_NC tmp_$CRESIS_FILE_NC
ncks -A usurf_$CRESIS_FILE_NC tmp_$CRESIS_FILE_NC
ncks -A topg_$CRESIS_FILE_NC tmp_$CRESIS_FILE_NC
ncatted -a units,Band1,o,c,"m" -a units,topg,o,c,"m" -a units,usurf,o,c,"m" tmp_$CRESIS_FILE_NC
nc2cdo.py --srs '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m' tmp_$CRESIS_FILE_NC
if [ [$NN == 1] ] ; then
  REMAP_EXTRAPOLATE=on cdo remapbil,$FL_FILE_NC tmp_$CRESIS_FILE_NC tmp2_$CRESIS_FILE_NC
else
  REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,$FL_FILE_NC tmp_$CRESIS_FILE_NC tmp2_$CRESIS_FILE_NC
fi
ncks -A -v x,y,mapping $FL_FILE_NC tmp2_$CRESIS_FILE_NC
ncatted -a grid_mapping,Band1,o,c,"mapping" tmp2_$CRESIS_FILE_NC

cdo -O setmisstoc,-9999. tmp2_$CRESIS_FILE_NC $CRESIS_FILE_NC
ncatted -a _FillValue,,d,, $CRESIS_FILE_NC
ncks -A -v x,y,mapping $FL_FILE_NC $CRESIS_FILE_NC
ncatted -a grid_mapping,Band1,o,c,"mapping" $CRESIS_FILE_NC

ncks -A -v Band1 $CRESIS_FILE_NC $FL_FILE_NC
ncatted -a _FillValue,,d,, $FL_FILE_NC
ncap2 -O -s "where(thk==-9999.) thk=Band1;" $FL_FILE_NC $FL_FILE_NC

ncrename -O -v Band1,thk $CRESIS_FILE_NC $CRESIS_FILE_NC

# GIMP DEM
GIMP=gimpdem_90m
GIMP_FILE_NC=${PROJECT}_gimp_${GS}m.nc
wget -nc ftp://ftp-bprc.mps.ohio-state.edu/downloads/gdg/gimpdem/$GIMP.tif.zip
unzip -o $GIMP.tif.zip
gdal_translate -projwin  $X_MIN $Y_MAX $X_MAX $Y_MIN $GIMP.tif ${PROJECT}_${GIMP}.tif
gdaldem hillshade -s 0.5 ${PROJECT}_${GIMP}.tif ${PROJECT}_${GIMP}_hillshade.tif
gdalwarp $WARPOPTIONS -of netCDF $GIMP.tif tmp_$GIMP_FILE_NC
ncrename -v Band1,usurf tmp_$GIMP_FILE_NC
nc2cdo.py --srs '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m' tmp_$GIMP_FILE_NC
if [ [$NN == 1] ] ; then
  REMAP_EXTRAPOLATE=on cdo remapbil,$FL_FILE_NC tmp_$GIMP_FILE_NC $GIMP_FILE_NC
else
  REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,$FL_FILE_NC tmp_$GIMP_FILE_NC $GIMP_FILE_NC
fi
ncks -A -v x,y,mapping $FL_FILE_NC $GIMP_FILE_NC
ncatted -a grid_mapping,usurf,o,c,"mapping"   $GIMP_FILE_NC

echo "Fetching University of Montana 1km data set ... "
UMT_FILE=Greenland1km.nc
UMT_FILE_NC=${PROJECT}_umt_${GS}m.nc
wget -nc http://websrv.cs.umt.edu/isis/images/a/ab/$UMT_FILE
nc2cdo.py $UMT_FILE
ncwa -O -a t $UMT_FILE tmp_$UMT_FILE
if [ [$NN == 1] ] ; then
  REMAP_EXTRAPOLATE=on cdo remapbil,$FL_FILE_NC -selvar,thk,topg tmp_$UMT_FILE $UMT_FILE_NC
else
  REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,$FL_FILE_NC tmp_$UMT_FILE $UMT_FILE_NC
fi
ncks -A -v x,y,mapping $FL_FILE_NC $UMT_FILE_NC
ncatted -a grid_mapping,thk,o,c,"mapping" -a grid_mapping,topg,o,c,"mapping" $UMT_FILE_NC
nc2cdo.py --srs '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m' $UMT_FILE_NC

# get SeaRISE file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
SR_FILE=Greenland_5km_v1.1.nc
SR_FILE_NC=${PROJECT}_searise_v1.1_${GS}m.nc

echo "Fetching SeaRISE master file ... "
wget -nc http://websrv.cs.umt.edu/isis/images/a/a5/$SR_FILE
echo "  ... done."
echo
nc2cdo.py $SR_FILE

# Regrid SeaRISE onto local grid
cdo selvar,smb,topg,thk $SR_FILE tmp_$SR_FILE
if [ [$NN == 1] ] ; then
  REMAP_EXTRAPOLATE=on cdo remapbil,$FL_FILE_NC tmp_$SR_FILE $SR_FILE_NC
else
  REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,$FL_FILE_NC tmp_$SR_FILE $SR_FILE_NC
fi
ncks -A -v x,y,mapping $FL_FILE_NC $SR_FILE_NC
ncatted -a grid_mapping,thk,o,c,"mapping" -a grid_mapping,topg,o,c,"mapping" -a grid_mapping,smb,o,c,"mapping" $SR_FILE_NC
nc2cdo.py --srs '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m' $SR_FILE_NC

# remap surface velocities, select area frist to speed
# things up a bit
VELIN_FILE=surf_vels.nc
cdo -O setmisstoc,0. $VELIN_FILE tmp_$VELIN_FILE
ncap2 -O -L 3 -4 -S removepoints.nco  tmp_$VELIN_FILE tmp_$VELIN_FILE
fill_missing.py -v magnitude,us,vs -e 1 -f tmp_$VELIN_FILE -o tmp2_$VELIN_FILE
ncks -A -v x,y,mapping $VELIN_FILE tmp2_$VELIN_FILE
MASK_FILE=g1km_0_CLRUN_mask.nc
python scripts/resample-mask.py -n $NN $MASK_FILE tmp2_$VELIN_FILE
ncap2 -O -s "where(mask==4) {magnitude=0.; us=0.; vs=0.;}" tmp2_$VELIN_FILE tmp2_$VELIN_FILE

VELOUT_FILE=${PROJECT}_surf_vels_${GS}m.nc

if [ [$NN == 1] ] ; then
  cdo remapbil,$FL_FILE_NC -sellonlatbox,$LON_MIN,$LON_MAX,$LAT_MIN,$LAT_MAX -selvar,us,vs,magnitude tmp2_$VELIN_FILE $VELOUT_FILE
else
  cdo -P $NN remapbil,$FL_FILE_NC  -sellonlatbox,$LON_MIN,$LON_MAX,$LAT_MIN,$LAT_MAX -selvar,us,vs,magnitude tmp2_$VELIN_FILE $VELOUT_FILE
fi
ncks -A -v x,y,mapping $FL_FILE_NC $VELOUT_FILE
ncatted -a grid_mapping,us,o,c,"mapping" -a grid_mapping,vs,o,c,"mapping" -a grid_mapping,magnitude,o,c,"mapping" $VELOUT_FILE
