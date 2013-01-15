#!/bin/bash

# You need to open a port first in a different shell:
# ssh -L5433:icebridge.sr.unh.edu:5432 bmachuser@icebridge.sr.unh.edu

# Kinda non-standard tools needed:
# CDO (https://code.zmaw.de/projects/cdo) e.g. from MacPorts
# Pyresample (https://code.google.com/p/pyresample/)
# nc2cdo.py (part of the PISM distribution)

set -x -e

# run ./prepare_jakobshavn.sh 1 if you havent CDO compiled with OpenMP
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
minlon=-51
maxlon=-42
minlat=68.25
maxlat=70

mod_val=1
mod_field='gid'

epsg_code_out=3413

filename=jak_basin

# Well this takes a while...
# But we leave it at that for now. In the longer run, we need a python 
# script that downloads the data and puts it into a more appriate format
# than ASCII.
# TODO:
# -9999 is the missing value. That should be dealt with in a smarter way.
# python scripts/general_query.py -table cresis_gr -fields "wgs84surf,wgs84bed" -epsg $epsg_code_out -and_clause "wgs84bed>-9999" -box $minlon $maxlon $minlat $maxlat  -mod_val $mod_val -mod_field gid > ${filename}.txt


# destination area in EPSG:3413 coordinates
xmin=-230000.0
xmax=80000.0
ymin=-2350000.0
ymax=-2200000.0
python scripts/preprocess.py -g $GS --bounds $xmin $xmax $ymin $ymax -n $NN ${filename}.txt tmp_${filename}.nc
# nc2cdo.py is from pism/util/
# it adds lat/lon, but also the 4 grid corners of each cell, needed for
# conservative remapping via CDO.
nc2cdo.py tmp_${filename}.nc
cdo setmisstoc,-9999. -selvar,thk tmp_${filename}.nc ${filename}.nc
ncatted -a _FillValue,,d,, ${filename}.nc
ncks -A -v thk -x tmp_${filename}.nc ${filename}.nc
WARPOPTIONS=" -overwrite -multi -r bilinear -te $xmin $ymin $xmax $ymax -tr $GS $GS -t_srs EPSG:$epsg_code_out"

# CReSIS data set
CRESIS=Jakobshavn_2006_2012_Composite
CRESISNC=cresis_thk.nc
wget -nc --no-check-certificate https://data.cresis.ku.edu/data/grids/$CRESIS.zip
unzip -o $CRESIS.zip
gdalwarp $WARPOPTIONS -of netCDF $CRESIS/grids/jakobshavn_2006_2012_composite_thickness.txt tmp_$CRESISNC

cdo -O setmisstoc,-9999. tmp_$CRESISNC $CRESISNC
ncatted -a _FillValue,,d,, $CRESISNC

ncks -A -v Band1 $CRESISNC ${filename}.nc
ncatted -a _FillValue,,d,, ${filename}.nc
ncap2 -O -s "where(thk==-9999.) thk=Band1;" ${filename}.nc ${filename}.nc
# GIMP DEM
GIMP=gimpdem_90m
wget -nc ftp://ftp-bprc.mps.ohio-state.edu/downloads/gdg/gimpdem/$GIMP.tif.zip
unzip -o $GIMP.tif.zip
gdalwarp $WARPOPTIONS -of netCDF $GIMP.tif $GIMP.nc
ncrename -v Band1,usurf $GIMP.nc
nc2cdo.py --srs '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m' $GIMP.nc
if [ [$NN == 1] ] ; then
  REMAP_EXTRAPOLATE=on cdo remapbil,${filename}.nc $GIMP.nc jak_$GIMP.nc
else
  REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,${filename}.nc $GIMP.nc jak_$GIMP.nc
fi
# TODO combine with SeaRISE file

# get SeaRISE file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
DATAVERSION=1.1
DATAURL=http://websrv.cs.umt.edu/isis/images/a/a5/
DATANAME=Greenland_5km_v$DATAVERSION.nc

echo "Fetching SeaRISE master file ... "
wget -nc ${DATAURL}${DATANAME}
echo "  ... done."
echo
nc2cdo.py $DATANAME

# Regrid SeaRISE onto local Jakobshavn grid

SEARISENAME=Greenland_5km_v${DATAVERSION}_small.nc
cdo selvar,smb,topg $DATANAME $SEARISENAME
OUTFILE=jak_input_v${DATAVERSION}.nc
if [ [$NN == 1] ] ; then
  REMAP_EXTRAPOLATE=on cdo remapbil,${filename}.nc $SEARISENAME $OUTFILE
else
  REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,${filename}.nc $SEARISENAME $OUTFILE
fi

ncks -A -v usurf jak_$GIMP.nc $OUTFILE
ncks -A -v Band1 $CRESISNC $OUTFILE
ncrename -v Band1,thk $OUTFILE

# CDO drops x,y and mapping information. Re-add it.
ncks -A -v x,y,mapping ${filename}.nc $OUTFILE
# and re-add grid_mapping attribute
ncatted -a grid_mapping,thk,o,c,"mapping" -a grid_mapping,topg,o,c,"mapping" -a grid_mapping,smb,o,c,"mapping" $OUTFILE

# remap surface velocities, select area frist to speed
# things up a bit
VELIN=surf_vels.nc
VELOUT=jak_surf_vels.nc
if [ [$NN == 1] ] ; then
  cdo remapbil,${filename}.nc -sellonlatbox,$minlon,$maxlon,$minlat,$maxlat $VELIN $VELOUT
else
  cdo -P $NN remapbil,${filename}.nc -sellonlatbox,$minlon,$maxlon,$minlat,$maxlat $VELIN $VELOUT
fi
# CDO drops x,y and mapping information. Re-add it.
ncks -A -v x,y,mapping ${filename}.nc $VELOUT
# and re-add grid_mapping attribute
ncatted -a grid_mapping,us,o,c,"mapping" -a grid_mapping,vs,o,c,"mapping" -a grid_mapping,magnitude,o,c,"mapping" $VELOUT
