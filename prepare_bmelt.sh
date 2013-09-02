
echo "Preparing basal melt"
FILE=g1km_0_CLRUN-bmelt.nc
nc2cdo.py $FILE
FILE_NC=${PROJECT}_bmelt_${GS}m.nc
if [ [$NN == 1] ] ; then
  REMAP_EXTRAPOLATE=on cdo remapbil,$FL_FILE_NC -selvar,thk,topg $FILE $FILE_NC
else
  REMAP_EXTRAPOLATE=on cdo -P $NN remapbil,$FL_FILE_NC $FILE $FILE_NC
fi
ncks -A -v x,y,mapping $FL_FILE_NC $FILE_NC
ncatted -a grid_mapping,bmelt,o,c,"mapping" $FILE_NC
nc2cdo.py --srs '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m' $FILE_NC
python scripts/scalar_within_poly.py -s 600 -v bmelt ~/data/data_sets/GreenlandFlightlines/submarine_melt.shp $FILE_NC

