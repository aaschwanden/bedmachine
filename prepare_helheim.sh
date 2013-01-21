#!/bin/bash

# You need to open a port first in a different shell:
# ssh -L5433:icebridge.sr.unh.edu:5432 bmachuser@icebridge.sr.unh.edu

# Kinda non-standard tools needed:
# CDO (https://code.zmaw.de/projects/cdo) e.g. from MacPorts
# Pyresample (https://code.google.com/p/pyresample/)
# nc2cdo.py (part of the PISM distribution)

set -x -e

# run ./prepare.sh 1 if you havent CDO compiled with OpenMP
NN=8  # default number of processors
if [ $# -gt 0 ] ; then
  NN="$1"
fi

# grid spacing
GS=500
if [ $# -gt 1 ] ; then
  GS="$2"
fi


# corners of lon,lat box for query
LON_MIN=-39.5
LON_MAX=-37.5
LAT_MIN=66.2
LAT_MAX=67.0

MOD_VAL=1
MOD_FIELD='gid'

EPSG=3413

PROJECT=helheim
PROJECTC=Helheim
YEARA=2008
YEARE=2012
CRESIS_YEARS=${YEARA}_${YEARE}

# destination area in EPSG:3413 coordinates
X_MIN=250000.0
X_MAX=330000.0
Y_MIN=-2590000.0
Y_MAX=-2510000.0

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# below here, everything should be general.

FL_FILE_TXT=${PROJECT}_flighlines.txt

# Well this takes a while...
# But we leave it at that for now. In the longer run, we need a python 
# script that downloads the data and puts it into a more appriate format
# than ASCII.
# TODO:
# -9999 is the missing value. That should be dealt with in a smarter way.

#python scripts/general_query.py -table cresis_gr -fields "wgs84surf,wgs84bed" -epsg $EPSG -and_clause "wgs84bed>-9999" -box $LON_MIN $LON_MAX $LAT_MIN $LAT_MAX  -mod_val $MOD_VAL -mod_field $MOD_FIELD > $FL_FILE_TXT

FL_FILE_NC=${PROJECT}_flighlines_${GS}m.nc
#python scripts/preprocess.py -g $GS --bounds $X_MIN $X_MAX $Y_MIN $Y_MAX -n $NN $FL_FILE_TXT tmp_$FL_FILE_NC

# nc2cdo.py is from pism/util/
# it adds lat/lon, but also the 4 grid corners of each cell, needed for
# conservative remapping via CDO.

# TODO figure out if we really need this.
cdo setmisstoc,-9999. -selvar,thk tmp_$FL_FILE_NC $FL_FILE_NC
ncatted -a _FillValue,,d,, $FL_FILE_NC
ncks -A -v thk -x tmp_$FL_FILE_NC $FL_FILE_NC
nc2cdo.py $FL_FILE_NC


WARPOPTIONS="-overwrite -multi -r bilinear -te $X_MIN $Y_MIN $X_MAX $Y_MAX -tr $GS $GS -t_srs EPSG:$EPSG"


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
VELOUT_FILE=${PROJECT}_surf_vels_${GS}m.nc
if [ [$NN == 1] ] ; then
  cdo remapbil,$FL_FILE_NC -sellonlatbox,$LON_MIN,$LON_MAX,$LAT_MIN,$LAT_MAX $VELIN_FILE $VELOUT_FILE
else
  cdo -P $NN remapbil,$FL_FILE_NC  -sellonlatbox,$LON_MIN,$LON_MAX,$LAT_MIN,$LAT_MAX $VELIN_FILE $VELOUT_FILE
fi
ncks -A -v x,y,mapping $FL_FILE_NC $VELOUT_FILE
ncatted -a grid_mapping,us,o,c,"mapping" -a grid_mapping,vs,o,c,"mapping" -a grid_mapping,magnitude,o,c,"mapping" $VELOUT_FILE
