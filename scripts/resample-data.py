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
parser.add_argument("--bounds", dest="bounds", nargs=4, type=float,
                    help="xmin xmax ymin ymax",
                    default=[-230000.0, 80000.0, -2350000.0, -2200000.0])
parser.add_argument("-g","--grid_spacing", dest="grid_spacing", type=float,
                  help='''target grid spacing in meters. Default=500m''', default=500)
parser.add_argument("-m","--missing_value", dest="miss_val", type=float,
                  help='''Missing value. Default=-9999.''', default=-9999.)
parser.add_argument("-n","--no_procs", dest="nprocs", type=int,
                  help='''No. of cores used for resamping.''', default=4)

options = parser.parse_args()
args = options.FILE
bounds = options.bounds
data_file = args[0]
if len(args) > 0:
    out_filename = args[1]
else:
    out_filename = "foo.nc"
grid_spacing = options.grid_spacing
miss_val = options.miss_val
nprocs = options.nprocs
fill_value = -2e-9
# FIXME: how to choose this value?
radius_of_influence = grid_spacing  # m

# read in data
# we could make this more flexible by only reading lon, lat and then use
# proj4 to convert to user-specified coordinate reference system.
x, y, lon, lat, thk = np.loadtxt(data_file, unpack=True, skiprows=5, delimiter=',')

xmin, xmax, ymin, ymax = bounds[0], bounds[1], bounds[2], bounds[3]
M = int((xmax - xmin) / grid_spacing)
N = int((ymax - ymin) / grid_spacing)
easting = np.linspace(xmin, xmax, M+1)
northing = np.linspace(ymin, ymax, N+1)
X, Y = np.meshgrid(easting, northing)

area_id = 'greenland'
area_name = 'Greenland'
proj_id = 'EPSG:3413'
proj4_str = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m'

area_extent = (xmin, ymin, xmax, ymax)
area_def = pr.utils.get_area_def(area_id, area_name, proj_id, proj4_str,
                              M+1, N+1, area_extent)


swath_def = pr.geometry.SwathDefinition(lons=lon, lats=lat)
result = pr.kd_tree.resample_nearest(swath_def, thk,
                                     area_def,
                                     radius_of_influence=radius_of_influence,
                                     epsilon=0.5,
                                     fill_value=fill_value,
                                     nprocs=nprocs)

# binary density (0: no data, 1: data)
rho = np.ones_like(result)
rho[result==fill_value] = 0

# Make a quick plot
fig_str = out_filename.split('.')[0]
fig_name = fig_str + '.png'
pr.plot.save_quicklook(fig_name, area_def, result, label='ice thickness')

# Write the data:
nc = CDF(out_filename, "w", format="NETCDF3_CLASSIC")
nc.createDimension("y", size=northing.shape[0])
nc.createDimension("x", size=easting.shape[0])
x = nc.createVariable("x", 'f', dimensions=("x",))
x.units = "m";
x.long_name = "easting"
x.standard_name = "projection_x_coordinate"

y = nc.createVariable("y", 'f', dimensions=("y",))
y.units = "m";
y.long_name = "northing"
y.standard_name = "projection_y_coordinate"

mapping_var = 'mapping'
mapping = nc.createVariable(mapping_var, 'b')
mapping.grid_mapping_name = "polar_stereographic"
mapping.latitude_of_projection_origin = 90.
mapping.straight_vertical_longitude_from_pole = -45.0
mapping.standard_parallel = 70.0
mapping.false_easting = 0.
mapping.false_northing = 0.
mapping.Northernmost_Northing = ymax
mapping.Southernmost_Northing = ymin
mapping.Easternmost_Easting = xmax
mapping.Westernmost_Easting = xmin
mapping.units = "m"

thk_var = nc.createVariable("thk", 'f', dimensions=("y", "x"),fill_value=fill_value)
thk_var.units = "m"
thk_var.long_name = "land ice thickness"
thk_var.standard_name = "land_ice_thickness"
thk_var.grid_mapping = mapping_var
thk_var[:] = np.flipud(result)

rho_var = nc.createVariable("rho", 'b', dimensions=("y", "x"),fill_value=fill_value)
rho_var.long_name = "data density integer mask"
rho_var.flag_values = "0b, 1b"
rho_var.flag_meanings = "no_data data"
rho_var.grid_mapping = mapping_var
rho_var[:] = np.flipud(rho)

x[:] = easting
y[:] = northing

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

