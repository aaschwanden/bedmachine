#!/usr/bin/env python

# Copyright (C) 2013 Andy Aschwanden

from argparse import ArgumentParser
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC
from osgeo import ogr

# Set up the option parser
parser = ArgumentParser()
parser.description = "A script for profile plots using pylab/matplotlib."
parser.add_argument("FILE", nargs=2)
parser.add_argument("-m", "--melt_rate",dest="melt_rate", type=float,
                  help="",default=900)

options = parser.parse_args()
args = options.FILE
melt_rate = options.melt_rate

import ogr
driver = ogr.GetDriverByName('ESRI Shapefile')
data_source = driver.Open(args[0], 0)
layer = data_source.GetLayer(0)
srs=layer.GetSpatialRef()
# Make sure we use lat/lon coordinates.
# Fixme: allow reprojection onto lat/lon if needed.
if not srs.IsGeographic():
    print('''Spatial Reference System in % s is not lat/lon. Exiting.'''
          % filename)
    import sys
    sys.exit(0)
feature = layer.GetFeature(0)

nc = NC(args[1], 'a')

bmelt = nc.variables['bmelt']
lat = nc.variables['lat']
lon = nc.variables['lon']

for m in range(0, bmelt.shape[1]):
    for n in range(0, bmelt.shape[2]):
        x = lon[m,n]
        y = lat[m,n]
        wkt = "POINT(%f %f)" % (x,y)
        point = ogr.CreateGeometryFromWkt(wkt)
        if feature.GetGeometryRef().Contains(point):
            bmelt[0,m,n] = melt_rate

nc.close()
