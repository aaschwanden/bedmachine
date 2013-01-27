#!/usr/bin/env python
# Copyright (C) 2013 Andy Aschwanden
#
# Example, using a 50m grid spacing and 8 cores:
# python scripts/preprocess.py -g 500 \
# --bounds -230000.0 80000.0 -2350000.0 -2200000.0 \
# -n 8 jakobshavn_flightlines.txt tmp_flightlines_500m.nc
# for the Jakobshavn basin, assuming the projection is EPSG:3413

import numpy as np
import pylab as plt
from argparse import ArgumentParser
from pyproj import Proj
import pyresample as pr
try:
    from netCDF3 import Dataset as CDF
except:
    from netCDF4 import Dataset as CDF

# Set up the option parser
parser = ArgumentParser()
parser.description = '''A script to preprocess ice thickness data.'''
parser.add_argument("FILE", nargs='*')
parser.add_argument("-n","--no_procs", dest="nprocs", type=int,
                  help='''No. of cores used for resamping.''', default=8)

options = parser.parse_args()
args = options.FILE
data_file = args[0]
if len(args) > 0:
    out_filename = args[1]
else:
    out_filename = "foo.nc"
nprocs = options.nprocs



nc = CDF(data_file, "r")
x_in = nc.variables['x'][:]
y_in = nc.variables['y'][:]
M_in = len(x_in)
N_in = len(y_in)
x_in_min = x_in[0]
x_in_max = x_in[-1]
y_in_min = y_in[0]
y_in_max = y_in[-1]
data = np.squeeze(nc.variables['mask'][:])
radius_of_influence = x_in[1] - x_in[0]  # m
nc.close()

area_id = 'greenland'
area_name = 'Greenland'
proj_id = 'Polar Stereo'
proj4_str = '+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m'

area_extent = (x_in_min, y_in_min, x_in_max, y_in_max)
area_in_def = pr.utils.get_area_def(area_id, area_name, proj_id, proj4_str,
                              M_in, N_in, area_extent)

# Write the data:
nc = CDF(out_filename, "a")
x_out = nc.variables['x'][:]
y_out = nc.variables['y'][:]
M_out = len(x_out)
N_out = len(y_out)
x_out_min = x_out[0]
x_out_max = x_out[-1]
y_out_min = y_out[0]
y_out_max = y_out[-1]


area_id = 'greenland'
area_name = 'Greenland'
proj_id = 'EPSG:3413'
proj4_str = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m'

area_extent = (x_out_min, y_out_min, x_out_max, y_out_max)
area_out_def = pr.utils.get_area_def(area_id, area_name, proj_id, proj4_str,
                              M_out, N_out, area_extent)


fill_value = 4  # ocean
area_con_in_nn = pr.image.ImageContainerNearest(data, area_in_def,
                                        radius_of_influence=radius_of_influence,
                                        fill_value=fill_value,
                                        nprocs=nprocs)

area_con_out_nn = area_con_in_nn.resample(area_out_def)
result = area_con_out_nn.image_data

pr.plot.save_quicklook('foo.png', area_out_def, result, label='mask')

try:
    var = nc.createVariable("mask", 'b', dimensions=("y", "x"))
except:
    var = nc.variables['mask']
#var[:] = np.flipud(result)
var[:] = result
# Save the projection information:
nc.projection = proj4_str

nc.Conventions = "CF-1.5"

# writing global attributes
import time
script_command = ' '.join([time.ctime(), ':', __file__.split('/')[-1],
                           ' '.join([str(l) for l in args])])
nc.history = script_command

print "writing to %s ...\n" % out_filename
print "run nc2cdo.py to add lat/lon variables" 
nc.close()

