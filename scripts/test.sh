#!bin/bash

# corners of lon,lat box for query
minlon=-52.5
maxlon=-36.5
minlat=67.45
maxlat=71.6

mod_val=1
mod_field='gid'

epsg_code_out=3413

ssh bmachuser@icebridge.sr.unh.edu python general/general_query.py -table cresis_gr -fields "wgs84surf,wgs84bed" -epsg $epsg_code_out -and_clause "wgs84surf\>-500" -box $minlon $maxlon $minlat $maxlat  -mod_val $mod_val -mod_field gid > jak.out

general_query.py -table cresis_gr -fields "wgs84surf,wgs84bed" -epsg 3413 -and_clause "wgs84surf>-500" -box -52.5 -36.5 67.45 71.1  -mod_val 25 -mod_field gid