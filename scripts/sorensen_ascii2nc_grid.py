#!/usr/bin/env python


from argparse import ArgumentParser
import numpy as np

from pyproj import Proj

try:
    import netCDF4 as netCDF
except:
    import netCDF3 as netCDF
CDF = netCDF.Dataset

from griddata import griddata

if __name__ == "__main__": 

    # Set up the option parser
    parser = ArgumentParser()
    parser.description = "Turns ascii into netCDF."
    parser.add_argument("FILE", nargs='*')
    parser.add_argument("--bounds", dest="bounds", nargs=4, type=float,
                        help="xmin xmax ymin ymax",
                        default=[-230000.0, 80000.0, -2350000.0, -2200000.0])
    parser.add_argument("--epsg", dest="epsg", type=int,
                        help="EPSG code",
                        default=3413)
    parser.add_argument("-m", "--mask_filename", dest="mask_file",
                      help="Name of the mask file.", default=None)
    parser.add_argument("-o", "--output_filename", dest="out_file",
                      help="Name of the output file.", default=None)
    parser.add_argument("-g", "--grid_spacing",dest="grid_spacing",
                      help="use X m grid spacing",
                      metavar="X",default=500)
    parser.add_argument("-s", "--stride", dest="stride",
                      help='''Use a stride for input data. Mostly used for debugging.''',
        default=1)

    options = parser.parse_args()
    args = options.FILE
    bounds = options.bounds
    grid_spacing = options.grid_spacing
    mask_file = options.mask_file
    out_file = options.out_file
    infile = args[0]
    fill_value = -2e9
    epsg = options.epsg
    stride = int(options.stride)

    # define output grid
    e0 = bounds[0]
    n0 = bounds[2] 
    e1 = bounds[1]
    n1 = bounds[3]

    de = dn =  int(grid_spacing)  # m
    M = int((e1 - e0) / de) + 1
    N = int((n1 - n0) / dn) + 1

    easting  = np.linspace(e0, e1, M)
    northing = np.linspace(n0, n1, N)
    ee, nn = np.meshgrid(easting,northing)

    proj = Proj(init="epsg:%s" % (epsg))

    lon, lat = proj(ee, nn, inverse=True)

    fromto = infile.split('.')[0].split('_')[0:2]
    year_start = fromto[0]
    year_end = fromto[1]
    years = int(year_end) - int(year_start)
    if not out_file:
        out_file = ("icesat_dhdt_%s_%s.nc" % (year_start, year_end))

    nc = CDF(out_file,'w')

    nc.createDimension("x", size=easting.shape[0])
    nc.createDimension("y", size=northing.shape[0])

    var = 'x'
    var_out = nc.createVariable(var, 'f', dimensions=("x"))
    var_out.axis = "X"
    var_out.long_name = "X-coordinate in Cartesian system"
    var_out.standard_name = "projection_x_coordinate"
    var_out.units = "meters"
    var_out[:] = easting

    var = 'y'
    var_out = nc.createVariable(var, 'f', dimensions=("y"))
    var_out.axis = "Y"
    var_out.long_name = "Y-coordinate in Cartesian system"
    var_out.standard_name = "projection_y_coordinate"
    var_out.units = "meters";
    var_out[:] = northing

    var = 'lon'
    var_out = nc.createVariable(var, 'f', dimensions=("y","x"))
    var_out.units = "degrees_east";
    var_out.valid_range = -180., 180.
    var_out.standard_name = "longitude"
    var_out[:] = lon

    var = 'lat'
    var_out = nc.createVariable(var, 'f', dimensions=("y","x"))
    var_out.units = "degrees_north";
    var_out.valid_range = -90., 90.
    var_out.standard_name = "latitude"
    var_out[:] = lat

    mapping = nc.createVariable("mapping",'c')
    mapping.ellipsoid = "WGS84"
    mapping.false_easting = 0.
    mapping.false_northing = 0.
    mapping.grid_mapping_name = "polar_stereographic"
    mapping.latitude_of_projection_origin = 90.
    mapping.standard_parallel = 71.
    mapping.straight_vertical_longitude_from_pole = -39.

    print("Reading file %s" % infile)
    lon_in, lat_in, dhdt = np.loadtxt(infile, unpack=True)

    x_in, y_in = proj(lon_in[::stride], lat_in[::stride])
    print("Interpolating using griddata...")
    Dhdt = griddata(x_in, y_in, dhdt[::stride], ee, nn, ext=0, nul=fill_value)
    print("...done.")

    if mask_file:
        nc_mask = CDF(mask_file, 'r')
        mask = np.squeeze(nc_mask.variables['ftt_mask'][:])
        if (mask.shape == Dhdt.shape):
            Dhdt[mask==0] = fill_value
        else:
            print("shape mismatch: %s and %s"
                  % (str(mask.shape), str(Dhdt.shape)))
        nc_mask.close()
    
    var = nc.createVariable("dhdt", 'f', dimensions=("y", "x"), fill_value=fill_value)
    var.units = "m year-1"
    var.mapping = "mapping"
    var.long_name = ("ice surface elevation change between %s and %s"
                     % (year_start, year_end)) 
    var.Source = '''Provided by L. Sandberg Sorensen.'''
    var[:] = Dhdt / years

    from time import asctime
    historystr = 'Created ' + asctime() + ' by Andy Aschwanden, University of Alaska Fairbanks, USA \n'
    nc.history = historystr
    nc.projection = proj.srs
    nc.Conventions = 'CF-1.5'

    nc.close()
