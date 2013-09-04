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
GS=250
if [ $# -gt 1 ] ; then
  GS="$2"
fi

EPSG=3413

YEAR=2008

PROJECT=jakobshavn
PROJECTC=Jakobshavn
CRESIS_YEARS=2006_2012

# destination area in EPSG:3413 coordinates

# large
#X_MIN=-230000.0
#X_MAX=80000.0
#Y_MIN=-2350000.0
#Y_MAX=-2200000.0

# medium
#X_MIN=-230000.0
#X_MAX=-100000.0
#Y_MIN=-2330000.0
#Y_MAX=-2205000.0

# small
X_MIN=-200000.0
X_MAX=-125000.0
Y_MIN=-2290000.0
Y_MAX=-2255000.0

FL_FILE_TXT=Jakobshavn_2007_2012_Composite_Flightlines_selected2.csv

FL_FILE_NC=${PROJECT}_flightlines_${YEAR}_${GS}m.nc
FL1985_FILE_NC=${PROJECT}_flightlines_1985_${GS}m.nc
python scripts/resample-cresis-data.py -g $GS --bounds $X_MIN $X_MAX $Y_MIN $Y_MAX \
    -c terminus_thickness.csv -n $NN $FL_FILE_TXT tmp_$FL_FILE_NC
nc2cdo.py --srs '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m' tmp_$FL_FILE_NC

cdo setmisstoc,-9999. -selvar,thk tmp_$FL_FILE_NC $FL_FILE_NC
ncatted -a _FillValue,,d,, $FL_FILE_NC
ncks -A -v thk -x tmp_$FL_FILE_NC $FL_FILE_NC
nc2cdo.py $FL_FILE_NC

MASK_FILE_NC=${PROJECT}_mask_${GS}m.nc
ncap2 -O -s "mask=thk*0;" $FL_FILE_NC $MASK_FILE_NC
python scripts/scalar_within_poly.py -s 1 -v mask ~/data/data_sets/GreenlandFlightlines/ice_thickness_zero.shp $MASK_FILE_NC
python scripts/scalar_within_poly.py -s 1 -v mask ~/data/data_sets/GreenlandFlightlines/ice_thickness_fjord.shp $MASK_FILE_NC
ncks -O -v thk,rho -x $MASK_FILE_NC $MASK_FILE_NC

 
#source prepare_velocities.sh
#source prepare_velocities_1985.sh
#source prepare_bmelt.sh


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
ncatted -a units,dhdt,o,c,"m year-1" $SPOT_FILE_NC

# CReSIS data set
CRESIS=${PROJECTC}_${CRESIS_YEARS}
CRESIS_FILE_ZIP=${CRESIS}_Composite.zip
CRESIS_FILE_NC=${PROJECT}_cresis_${YEAR}_${GS}m.nc
CRESIS1985_FILE_NC=${PROJECT}_cresis_1985_${GS}m.nc
#wget -nc --no-check-certificate https://data.cresis.ku.edu/data/grids/$CRESIS_FILE_ZIP
#unzip -o $CRESIS_FILE_ZIP
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
nccopy $CRESIS_FILE_NC $CRESIS1985_FILE_NC


ncks -A -v Band1 $CRESIS_FILE_NC $FL_FILE_NC
ncatted -a _FillValue,,d,, $FL_FILE_NC
ncap2 -O -s "where(thk==-9999.) thk=Band1;" $FL_FILE_NC $FL_FILE_NC

ncrename -O -v Band1,thk $CRESIS_FILE_NC $CRESIS_FILE_NC


IN_DEM=DEM_5.5_july_24_85.nc
USURF1985_FILE_NC=${PROJECT}_usurf_1985_${GS}m.nc
if [ [$NN == 1] ] ; then
  REMAP_EXTRAPOLATE=on cdo remapbil,$FL_FILE_NC $IN_DEM tmp_$USURF1985_FILE_NC
else
  REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,$FL_FILE_NC $IN_DEM tmp_$USURF1985_FILE_NC
fi

ncks -A -v x,y,mapping $FL_FILE_NC tmp_$USURF1985_FILE_NC
ncatted -a grid_mapping,usurf,o,c,"mapping" -a units,usurf,o,c,"m"  -a long_name,usurf,o,c,"ice upper surface elevation" -a standard_name,usurf,o,c,"surface_altitude" tmp_$USURF1985_FILE_NC
fill_missing.py -v usurf -e 1.5 -f tmp_$USURF1985_FILE_NC -o $USURF1985_FILE_NC

USURF2008_FILE_NC=${PROJECT}_usurf_${YEAR}_${GS}m.nc
ncks -O $CRESIS_FILE_NC $USURF2008_FILE_NC
nc2cdo.py $USURF2008_FILE_NC

USURFDIFF_FILE_NC=${PROJECT}_usurf_2008-1985_${GS}m.nc
cdo sub -selvar,usurf $USURF2008_FILE_NC -selvar,usurf $USURF1985_FILE_NC $USURFDIFF_FILE_NC
ncks -A -v x,y,mapping $FL_FILE_NC $USURFDIFF_FILE_NC
ncatted -a grid_mapping,usurf,o,c,"mapping" $USURFDIFF_FILE_NC
nc2cdo.py $USURFDIFF_FILE_NC

cdo sub -selvar,usurf $CRESIS1985_FILE_NC -selvar,usurf $USURFDIFF_FILE_NC tmp_$CRESIS1985_FILE_NC

ncrename -v usurf,thk $USURFDIFF_FILE_NC
cdo sub -selvar,thk $CRESIS_FILE_NC -selvar,thk $USURFDIFF_FILE_NC tmp2_$CRESIS1985_FILE_NC
ncks -A -v thk tmp2_$CRESIS1985_FILE_NC $CRESIS1985_FILE_NC
ncks -A -v usurf tmp_$CRESIS1985_FILE_NC $CRESIS1985_FILE_NC
ncks -A -v x,y,mapping $FL_FILE_NC $CRESIS1985_FILE_NC
nc2cdo.py $CRESIS1985_FILE_NC

cdo sub -selvar,thk $FL_FILE_NC $USURFDIFF_FILE_NC tmp_${FL1985_FILE_NC}
ncks -A -v x,y,mapping $FL_FILE_NC tmp_${FL1985_FILE_NC}
ncatted -a grid_mapping,thk,o,c,"mapping" tmp_${FL1985_FILE_NC}
nccopy $FL_FILE_NC $FL1985_FILE_NC
ncks -A -v thk tmp_${FL1985_FILE_NC} $FL1985_FILE_NC
ncatted -a _FillValue,thk,o,f,-2e9 $FL1985_FILE_NC

python scripts/scalar_within_poly.py -s 0 -v thk ~/data/data_sets/GreenlandFlightlines/ice_thickness_zero.shp $CRESIS_FILE_NC
python scripts/scalar_within_poly.py -s 750 -v thk ~/data/data_sets/GreenlandFlightlines/ice_thickness_fjord.shp $CRESIS_FILE_NC

#source prepare_additional.sh
