#!/usr/bin/env python
# Copyright (C) 2013 Andy Aschwanden
#
# Example:
# python preprocess.py --bounds -241000 207000 -2384000 -2150000 jak_basin.txt
# for the Jakobshavn basin, assuming the projection is EPSG:3413

import numpy as np
import pylab as plt
from argparse import ArgumentParser
from griddata import griddata

# Set up the option parser
parser = ArgumentParser()
parser.description = "A script to plot a variable in a netCDF file over a GeoTiff. Uses GDAL python bindings, Proj4, and Basemap. Script is fine-tuned for whole Greenland plots, but can be adapted for other needs."
parser.add_argument("FILE", nargs='*')
parser.add_argument("--bounds", dest="bounds", nargs=4, type=float,
                    help="xmin xmax ymin ymax",
                    default=[-240000.0, 200000.0, -2360000.0, -2160000.0])
parser.add_argument("-g","--grid_spacing", dest="grid_spacing", type=float,
                  help='''target grid spacing in meters. Default=500m''', default=500)
parser.add_argument("-m","--missing_value", dest="miss_val", type=float,
                  help='''Missing value. Default=-9999.''', default=-9999.)

options = parser.parse_args()
args = options.FILE
bounds = options.bounds
data_file = args[0]
grid_spacing = options.grid_spacing
miss_val = options.miss_val


# read in data
# we could make this more flexible by only reading lon, lat and then use
# proj4 to convert to user-specified coordinate reference system.
x, y, lon, lat, wgs84surf, wgs84bed = np.loadtxt(data_file, unpack=True, skiprows=5, delimiter=',')

# find missing values
# FIXME: not very smart, there could be missing values in other fields too
idx = np.nonzero(wgs84bed==miss_val)

x = np.delete(x, idx)
y = np.delete(y, idx)
lat = np.delete(lat, idx)
lon = np.delete(lon, idx)
wgs84surf = np.delete(wgs84surf, idx)
wgs84bed = np.delete(wgs84bed, idx)

# also, we should do some quality check here
wgs84thk = (wgs84surf - wgs84bed)

xmin, xmax, ymin, ymax = bounds[0], bounds[1], bounds[2], bounds[3]
M = int((xmax - xmin) / grid_spacing)
N = int((ymax - ymin) / grid_spacing)
easting = np.linspace(xmin, xmax, M+1)
northing = np.linspace(ymin, ymax, N+1)
X, Y = np.meshgrid(easting, northing)

# just for debugging we define a skip
skip = 1

# Calculate gridded ice thickness using natural neighbor interpolation
# (http://code.google.com/p/griddata-python/)
Thk = griddata(x[::skip], y[::skip], wgs84thk[::skip], X, Y)

