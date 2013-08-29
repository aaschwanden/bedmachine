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

