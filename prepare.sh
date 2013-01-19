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
MIN_LON=-51
MAX_LON=-42
MIN_LAT=68.25
MAX_LAT=70

MOD_VAL=1
MOD_FIELD='gid'

EPSG=3413

PROJECT=jakobshavn



# destination area in EPSG:3413 coordinates
X_MIN=-230000.0
X_MAX=80000.0
Y_MIN=-2350000.0
Y_MAX=-2200000.0

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# below here, everything should be general.

FLFILE_TXT=${PROJECT}_flighlines.txt

# Well this takes a while...
# But we leave it at that for now. In the longer run, we need a python 
# script that downloads the data and puts it into a more appriate format
# than ASCII.
# TODO:
# -9999 is the missing value. That should be dealt with in a smarter way.
python scripts/general_query.py -table cresis_gr -fields "wgs84surf,wgs84bed" -epsg $EPSG -and_clause "wgs84bed>-9999" -box $MIN_LON $MAX_LON $MIN_LAT $MAX_LAT  -mod_val $MOD_VAL -mod_field MOD_FIELD > $FLFILE_TXT

FLFILE_NC=${PROJECT}_flighlines_${GS}.nc
python scripts/preprocess.py -g $GS --bounds $X_MIN $X_MAX $Y_MIN $Y_MAX -n $NN $FLFILE_TXT tmp_$FLFILE_NC
# nc2cdo.py is from pism/util/
# it adds lat/lon, but also the 4 grid corners of each cell, needed for
# conservative remapping via CDO.

# TODO figure out if we really need this.
cdo setmisstoc,-9999. -selvar,thk tmp_$FLFILE_NC $FLFILE_NC
ncatted -a _FillValue,,d,, $FLFILE_NC
ncks -A -v thk -x tmp_$FLFILE_NC $FLFILE_NC
nc2cdo.py $FLFILE_NC


WARPOPTIONS="-overwrite -multi -r bilinear -te $X_MIN $Y_MIN $X_MAX $Y_MAX -tr $GS $GS -t_srs EPSG:$EPSG"

exit
# CReSIS data set
CRESIS=${PROJECT}_${YEARS}_Composite
CRESISNC=cresis_thk_${GS}m.nc
wget -nc --no-check-certificate https://data.cresis.ku.edu/data/grids/$CRESIS.zip
unzip -o $CRESIS.zip
gdalwarp $WARPOPTIONS -of netCDF $CRESIS/grids/jakobshavn_2006_2012_composite_thickness.txt tmp_$CRESISNC
nc2cdo.py --srs '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m' tmp_$CRESISNC
if [ [$NN == 1] ] ; then
  REMAP_EXTRAPOLATE=on cdo remapbil,${filename}_${GS}m.nc tmp_$CRESISNC tmp2_$CRESISNC
else
  REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,${filename}_${GS}m.nc tmp_$CRESISNC tmp2_$CRESISNC
fi

cdo -O setmisstoc,-9999. tmp2_$CRESISNC $CRESISNC
ncatted -a _FillValue,,d,, $CRESISNC

ncks -A -v Band1 $CRESISNC ${filename}_${GS}m.nc
ncatted -a _FillValue,,d,, ${filename}_${GS}m.nc
ncap2 -O -s "where(thk==-9999.) thk=Band1;" ${filename}_${GS}m.nc ${filename}_${GS}m.nc
# GIMP DEM
GIMP=gimpdem_90m
wget -nc ftp://ftp-bprc.mps.ohio-state.edu/downloads/gdg/gimpdem/$GIMP.tif.zip
unzip -o $GIMP.tif.zip
gdalwarp $WARPOPTIONS -of netCDF $GIMP.tif $GIMP.nc
ncrename -v Band1,usurf $GIMP.nc
nc2cdo.py --srs '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m' $GIMP.nc
if [ [$NN == 1] ] ; then
  REMAP_EXTRAPOLATE=on cdo remapbil,${filename}_${GS}m.nc $GIMP.nc jak_${GIMP}_${GS}m.nc
else
  REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,${filename}_${GS}m.nc $GIMP.nc jak_${GIMP}_${GS}m.nc
fi

echo "Fetching University of Montana 1km data set ... "
wget -nc http://websrv.cs.umt.edu/isis/images/a/ab/Greenland1km.nc

# get SeaRISE file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
DATAVERSION=1.1
DATAURL=http://websrv.cs.umt.edu/isis/images/a/a5/
DATANAME=Greenland_5km_v$DATAVERSION.nc


echo "Fetching SeaRISE master file ... "
wget -nc ${DATAURL}${DATANAME}
echo "  ... done."
echo
nc2cdo.py $DATANAME

# Regrid SeaRISE onto local grid

SEARISENAME=Greenland_5km_v${DATAVERSION}_small.nc
cdo selvar,smb,topg $DATANAME $SEARISENAME
OUTFILE=jak_input_v${DATAVERSION}_${GS}m.nc
if [ [$NN == 1] ] ; then
  REMAP_EXTRAPOLATE=on cdo remapbil,${filename}_${GS}m.nc $SEARISENAME $OUTFILE
else
  REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,${filename}_${GS}m.nc $SEARISENAME $OUTFILE
fi

ncks -A -v usurf jak_${GIMP}_${GS}m.nc $OUTFILE
ncks -A -v Band1 $CRESISNC $OUTFILE
ncrename -v Band1,thk $OUTFILE

# CDO drops x,y and mapping information. Re-add it.
ncks -A -v x,y,mapping ${filename}_${GS}m.nc $OUTFILE
# and re-add grid_mapping attribute
ncatted -a grid_mapping,thk,o,c,"mapping" -a grid_mapping,topg,o,c,"mapping" -a grid_mapping,smb,o,c,"mapping" $OUTFILE

# remap surface velocities, select area frist to speed
# things up a bit
VELIN=surf_vels.nc
VELOUT=jak_surf_vels_${GS}m.nc
if [ [$NN == 1] ] ; then
  cdo remapbil,${filename}_${GS}m.nc -sellonlatbox,$minlon,$maxlon,$minlat,$maxlat $VELIN $VELOUT
else
  cdo -P $NN remapbil,${filename}_${GS}m.nc -sellonlatbox,$minlon,$maxlon,$minlat,$maxlat $VELIN $VELOUT
fi
# CDO drops x,y and mapping information. Re-add it.
ncks -A -v x,y,mapping ${filename}_${GS}m.nc $VELOUT
# and re-add grid_mapping attribute
ncatted -a grid_mapping,us,o,c,"mapping" -a grid_mapping,vs,o,c,"mapping" -a grid_mapping,magnitude,o,c,"mapping" $VELOUT
