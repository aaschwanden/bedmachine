#!/usr/bin/env python
import time
from shutil import copy
import numpy as np
from argparse import ArgumentParser

try:
    from netCDF3 import Dataset as CDF
except:
    from netCDF4 import Dataset as CDF

def stencil_4(i, j):
    '''
    Return a 4-star stencil for * (i,j)
            (i,j+1)
    (i-1,j)    *    (i+1,j)
            (i,j-1)
    '''
    return [i+1,i,i-1,i], [j,j+1,j,j-1]

# Set up the argument parser
parser = ArgumentParser()
parser.description ='''
A script to interpolate velocities in areas where the magnitude of errors
is above a given threshold. Direction and magnitude are averaged.
'''
parser.add_argument("FILE", nargs='*')
parser.add_argument("-d", "--debug", dest="debug", action="store_true",
                  help="Makes plots for debugging", default=False)
parser.add_argument("-e","--error_threshold", dest="error_threshold", type=float,
                    help='''Magnitude of error threshold above which interpolation is done''',
                    default=100.)
parser.add_argument("-o","--output_filename", dest="output_filename",
                    help='''Name of output file''',
                    default='foo.nc')

options = parser.parse_args()
args = options.FILE
debug = options.debug
error_threshold = options.error_threshold
output_filename = options.output_filename

input_filename = args[0]
try:
    copy(input_filename, output_filename)
except:
    print("Could not copy file %s to %s" % (input_filename, output_filename))

try:
    nc = CDF(output_filename, 'a')
except:
    print("Could not open file %s" % output_filename)


print(("    - reading variable %s from file %s" % ("us", output_filename)))
try:
    us = nc.variables["us"]
    us_missing = us._FillValue
except:
    print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
          % ("us", output_filename)))
    import sys
    sys.exit(1)

print(("    - reading variable %s from file %s" % ("vs", output_filename)))
try:
    vs = nc.variables["vs"]
    vs_missing = vs._FillValue
except:
    print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
          % ("vs", output_filename)))
    import sys
    sys.exit(1)

print(("    - reading variable %s from file %s" % ("ue", output_filename)))
try:
    ue = nc.variables["ue"][:]
except:
    print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
          % ("ue", output_filename)))
    import sys
    sys.exit(1)

print(("    - reading variable %s from file %s" % ("ve", output_filename)))
try:
    ve = nc.variables["ve"][:]
except:
    print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
          % ("ve", output_filename)))
    import sys
    sys.exit(1)

print(("    - reading variable %s from file %s" % ("magnitude", output_filename)))
try:
    magnitude = nc.variables["magnitude"]
    mag_missing = magnitude._FillValue
except:
    print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
          % ("magnitude", output_filename)))
    import sys
    sys.exit(1)

u = us[:]
v = vs[:]
Umag = magnitude[:]
error_mag = ue

if debug:
    import pylab as plt
    from matplotlib import colors
    norm = colors.LogNorm(vmin=np.min(error_mag), vmax=np.max(error_mag))
    plt.imshow(error_mag, norm=norm)

e_i, e_j = np.nonzero(error_mag > error_threshold)
pthresh = 3  # min 3 values required
# basis vector to calculate angle
e1 = [1, 0]

for i, j in zip(e_i, e_j):
    ii, jj = stencil_4(i,j)
    val = magnitude[ii,jj]
    S = np.sum(val != mag_missing)
    if S >= pthresh:
        mag = np.sum(val[val != mag_missing]) / S
        U = [u[i,j],v[i,j]]
        cosphi = np.dot(U, e1) / (np.linalg.norm(U,2) * np.linalg.norm(e1,2))
        if U[1] > 0:
            phi = np.arccos(cosphi)
        else:
            phi = 2*np.pi - np.arccos(cosphi)
        sinphi = np.sin(phi)
        Umag[i,j] = mag
        u[i,j] = mag * cosphi
        v[i,j] = mag * sinphi
    else:
        Umag[i,j] = mag_missing
        u[i,j] = us_missing
        v[i,j] = vs_missing

nc.close()

if debug:
    plt.show()
