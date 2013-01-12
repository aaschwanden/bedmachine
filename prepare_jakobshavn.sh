#!bin/bash

# to connect via QGIS postigs, you need to open a port first:
# ssh -L5433:icebridge.sr.unh.edu:5432 bmachuser@icebridge.sr.unh.edu
# and then conncect to "localhost" with the above port

# corners of lon,lat box for query
minlon=-51
maxlon=-40
minlat=68.25
maxlat=70.25

mod_val=1
mod_field='gid'

epsg_code_out=3413

ssh bmachuser@icebridge.sr.unh.edu python general/general_query.py -table cresis_gr -fields "wgs84surf,wgs84bed" -epsg $epsg_code_out -and_clause "wgs84surf\>-500" -box $minlon $maxlon $minlat $maxlat  -mod_val $mod_val -mod_field gid > jak_basin.txt
