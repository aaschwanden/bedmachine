#!/usr/bin/env python

# Copyright (C) 2012-2013 Jesse Johnson
#
# Modifications by Andy Aschwanden

import numpy as np
try:
    from netCDF3 import Dataset as CDF
except:
    from netCDF4 import Dataset as CDF
from dolfin import *

from argparse import ArgumentParser


class DataOutput:
    def __init__(self,directory):
        if directory[-1]=='/':
            self.directory = directory
        else:
            self.directory = directory+'/'
            
    def write_dictionary_of_files(self,d,extension='pvd'):
        """ Looking for a dictionary d of data to save. The keys are the file 
        names, and the values are the data fields to be stored. Also takes an
        optional extension to determine if it is pvd or xml output."""
        for filename in d:
            file_handle = File(self.directory+filename+'.'+extension)
            file_handle<<d[filename]
    def write_one_file(self,name,data,extension='pvd'):
        file_handle = File(self.directory+name+'.'+extension)
        file_handle<<data



def generate_expression_from_gridded_data(x, y, var, method='bil'):
    from scipy.interpolate import RectBivariateSpline
    from scipy.interpolate import NearestNDInterpolator
    if (method=="bil"):
      interpolant = RectBivariateSpline(x, y, var)
    elif (method=="nearest"):
      X,Y = np.meshgrid(x,y)
      coords = np.vstack((X.flatten(),Y.flatten())).T
      interpolant = NearestNDInterpolator(coords, var.flatten())
    else:
      print("method not recongnized: %s" % method)
    class newExpression(Expression):
        def __init_(self):
          pass
        def eval(self,values,x):
            values[0] = interpolant(x[0], x[1])

    return newExpression()


def get_dims(nc):
    '''
    Gets dimensions from netcdf instance

    Parameters:
    -----------
    nc: netCDF instance

    Returns:
    --------
    xdim, ydim, zdim, tdim: dimensions
    '''
        
    ## a list of possible x-dimensions names
    xdims = ['x','x1']
    ## a list of possible y-dimensions names
    ydims = ['y','y1']
    ## a list of possible z-dimensions names
    zdims= ['z', 'z1']
    ## a list of possible time-dimensions names
    tdims= ['t', 'time']

    ## assign x dimension
    for dim in xdims:
        if dim in list(nc.dimensions.keys()):
            xdim = dim
    ## assign y dimension
    for dim in ydims:
        if dim in list(nc.dimensions.keys()):
            ydim = dim
    ## assign y dimension
    for dim in zdims:
        if dim in list(nc.dimensions.keys()):
            zdim = dim
        else:
            zdim = 'z'
    ## assign y dimension
    for dim in tdims:
        if dim in list(nc.dimensions.keys()):
            tdim = dim
        else:
            tdim = 'time'
    return xdim, ydim, zdim, tdim


def permute(variable, output_order=('time', 'z', 'zb', 'y', 'x')):
    '''
    Permute dimensions of a NetCDF variable to match the output
    storage order.

    Parameters
    ----------
    variable : a netcdf variable
               e.g. thk = nc.variables['thk']
    output_order: dimension tuple (optional)
                  default ordering is ('time', 'z', 'zb', 'y', 'x')

    Returns
    -------
    var_perm : array_like
    '''

    input_dimensions = variable.dimensions

    # filter out irrelevant dimensions
    dimensions = filter(lambda(x): x in input_dimensions,
                        output_order)

    # create the mapping
    mapping = map(lambda(x): dimensions.index(x),
                  input_dimensions)

    if mapping:
        return np.transpose(variable[:], mapping)
    else:
        return variable[:]  # so that it does not break processing "mapping"


# Set up the option parser
parser = ArgumentParser()
parser.description = '''A script to preprocess ice thickness data.'''
parser.add_argument("FILE", nargs='*')
parser.add_argument("-s","--scale_factor", dest="scale_factor",
                    type=float,
                    help='''Scale the element size Values<1 mean a coarser mesh. Default=1''',
                    default=1.)
parser.add_argument("-a","--alpha", dest="alpha", type=float,
                    help='''Regularization parameter''',
                    default=2.0)
parser.add_argument("-g","--gamma", dest="gamma", type=float,
                    help='''Misfit penalty''',
                    default=2.0)
parser.add_argument("-i","--in_file", dest="data_filename",
                    help='''File containing target thickness and density''',
                    default='jak_basin.nc')
parser.add_argument("-v","--velocity_file", dest="vel_filename",
                    help='''File containing x,y components of surface velocity''',
                    default='jak_surface_vels.nc')
parser.add_argument("-b","--bc_file", dest="bc_filename",
                    help='''File containing ice surface, and other boundary conditions.''',
                    default='jak_input_v1.1..nc')
options = parser.parse_args()
scale_factor = options.scale_factor
parameters['allow_extrapolation'] = True
data_filename = options.data_filename
vel_filename = options.vel_filename
bc_filename = options.bc_filename
# Misfit penalty.  Penalized differences between the calculated and observed thickness.
gamma = options.gamma
# Regularization parameter (penalty on the gradient of the solution)
alpha = options.alpha

# minimum ice thickness
thk_min = 10
output_order = ("x", "y")


nc_data = CDF(data_filename, 'r')
xdim, ydim, zdim, tdim = get_dims(nc_data)
x = nc_data.variables[xdim][:]
y = nc_data.variables[ydim][:]
input_dimensions = nc_data.variables["thk"].dimensions
rho = np.squeeze(permute(nc_data.variables["rho"], output_order=output_order))
Hfl = np.squeeze(permute(nc_data.variables["thk"], output_order=output_order))
M = len(x)
xmin = x[0]
xmax = x[-1]
N = len(y)
ymin = y[0]
ymax = y[-1]
nc_data.close()

nc_vel = CDF(vel_filename, 'r')
xdim, ydim, zdim, tdim = get_dims(nc_vel)
uvel = np.squeeze(permute(nc_vel.variables["us"], output_order=output_order))
uvel_fill = nc_vel.variables["us"]._FillValue
uvel[uvel.data==uvel_fill] = 0.
vvel = np.squeeze(permute(nc_vel.variables["vs"], output_order=output_order))
vvel_fill = nc_vel.variables["vs"]._FillValue
vvel[vvel.data==vvel_fill] = 0.
nc_vel.close()

output_order = ("time", "x", "y")
nc_bc = CDF(bc_filename, 'r')
xdim, ydim, zdim, tdim = get_dims(nc_bc)
Hobs = np.squeeze(permute(nc_bc.variables["thk"], output_order=output_order))
Hobs[Hobs<thk_min] = thk_min
S = np.squeeze(permute(nc_bc.variables["usurf"], output_order=output_order))
smb = np.squeeze(permute(nc_bc.variables["smb"], output_order=output_order))
nc_bc.close()

set_log_level(PROGRESS)

extend = 50
MM = int(np.ceil(M * scale_factor))
NN = int(np.ceil(N * scale_factor))
x_scaled = np.linspace(xmin - extend, xmax + extend, MM)
y_scaled = np.linspace(ymin - extend, ymax + extend, NN)
mesh = RectangleMesh(np.float(xmin) - extend, np.float(ymin - extend),
                        np.float(xmax + extend) ,np.float(ymax + extend),
                        MM,
                        NN
                        )

func_space =  FunctionSpace(mesh, "CG", 1)
func_space_dg =  FunctionSpace(mesh, "DG", 1)

Hin = project(generate_expression_from_gridded_data(x, y, Hobs), func_space)
H0 = project(generate_expression_from_gridded_data(x, y, Hfl), func_space)
S = project(generate_expression_from_gridded_data(x, y, S), func_space)
adot = project(generate_expression_from_gridded_data(x, y, smb), func_space)
rho_d = project(generate_expression_from_gridded_data(x, y, rho, method="nearest"),
                func_space_dg)
u_o = project(generate_expression_from_gridded_data(x, y, uvel), func_space)
v_o = project(generate_expression_from_gridded_data(x, y, vvel), func_space)

# Velocity norm
U = as_vector([u_o,v_o])  # unsmoothed
Unorm = project(sqrt(dot(U,U)) + 1e-10)

# Steepest descents velocities (if needed)
u_s = project(-S.dx(0) * Unorm)
v_s = project(-S.dx(1) * Unorm)

# Ignore slow regions of ice.
utol = 5.0

# Is this correct??
def inside(x, on_boundary):
  return (Unorm(x[0],x[1]) < utol) or (on_boundary and \
                       (x[0] < DOLFIN_EPS or x[1] < DOLFIN_EPS or \
                       (x[0] > 0.5 - DOLFIN_EPS and x[1] > 0.5 - DOLFIN_EPS)))

dbc = DirichletBC(func_space, Hin, inside)

# Solution and Trial function
H = Function(func_space)
dH = TrialFunction(func_space)

# Test Function
phi = TestFunction(func_space)

# Objective function
I = (div(U*H) - adot)**2*dx \
    + gamma*rho_d*0.5*(H-H0)**2*dx \
    + alpha*(H.dx(0)**2 + H.dx(1)**2)*dx

delta_I = derivative(I, H, phi)

J = derivative(delta_I, H, dH)

solve(delta_I==0, H, dbc, J=J)

prefix, suffix = data_filename.split('.')
# Dolfin file handler doesn't like unicode
prefix = str(prefix)
gamma_str = '_'.join(['gamma', str(gamma)])
alpha_str = '_'.join(['alpha', str(alpha)])

do = DataOutput('./')
data_out = {'_'.join([prefix, alpha_str, gamma_str, 'mcb_bed']) : project(S-H),
            '_'.join([prefix, alpha_str, gamma_str, 'cresis_bed']) : project(S-Hin),
            '_'.join([prefix, alpha_str, gamma_str, 'cresis_flux_div_obs']) : project(div(U*Hin)),
            '_'.join([prefix, alpha_str, gamma_str, 'mcb_flux_div']) : project(div(U*H)),
            '_'.join([prefix, alpha_str, gamma_str, 'U']) : Unorm,
            '_'.join([prefix, alpha_str, gamma_str, 'rho_d']) : rho_d,
            '_'.join([prefix, alpha_str, gamma_str, 'adot']) : adot,
            '_'.join([prefix, alpha_str, gamma_str, 'S']) : S,
            '_'.join([prefix, alpha_str, gamma_str, 'H0']) : H0,
            '_'.join([prefix, alpha_str, gamma_str, 'thk']) : H,
            '_'.join([prefix, alpha_str, gamma_str, 'delta_H'])  : project(H-H0),
            '_'.join([prefix, alpha_str, gamma_str, 'delta_H_obs'])  : project(H-Hin)
            }
do.write_dictionary_of_files(data_out)


def evaluate_regular_grid(f, x, y):
    fa = np.zeros((y.size, x.size))
    parameters['allow_extrapolation']=True
    for i, xs in enumerate(x):
        for j, ys in enumerate(y):
            try:
                fa[j,i] = f(xs,ys)
            except:
                print "Point not in mesh, skipping."
                pass
    return fa

from shutil import copy, move
from tempfile import mkstemp
from os import close

input_filename = data_filename
output_filename = '_'.join([prefix, alpha_str, gamma_str]) + '.nc'

## mesh_out = RectangleMesh(np.float(xmin), np.float(ymin),
##                          np.float(xmax), np.float(ymax),
##                          M, N)

## func_space_out =  FunctionSpace(mesh_out, "CG", 1)
## func_space_dg_out =  FunctionSpace(mesh_out, "DG", 1)

print "Creating the temporary file..."
try:
    (handle, tmp_filename) = mkstemp()
    close(handle) # mkstemp returns a file handle (which we don't need)
    copy(input_filename, tmp_filename)
except IOError:
    print "ERROR: Can't create %s, Exiting..." % tmp_filename

try:
    nc = CDF(tmp_filename, 'a')
except Exception, message:
   print message
   print "Note: %s was not modified." % output_filename
   exit(-1)


def create_variable(var_name, f, long_name=None, standard_name=None, units=None, fill_value=None):
    var = nc.createVariable(var_name, np.double, input_dimensions)
    if long_name is not None:
        var.long_name = long_name
    if standard_name is not None:
        var.standard_name = standard_name
    if units is not None:
        var.units = units
    var.grid_mapping = "mapping"
    if fill_value is not None:
        var._FillValue = fill_value
    var[:] = evaluate_regular_grid(f, x, y)



create_variable("thk_mcb", project(H), long_name="land ice thickness from mass conservation", units="m")
create_variable("delta_thk", project(H-H0), long_name="difference", units="m")
create_variable("topg_mcb", project(S-H), long_name="bedrock surface elevation from mass conservation", standard_name="bedrock_altitude", units="m")
create_variable("divHU_mcb", project(div(H*U)), long_name="flux divergence using mass conservation", units="m year-1")
create_variable("divHU_in", project(div(Hin*U)), long_name="flux divergence observed", units="m year-1")

## H_out = interpolate(H, func_space_out)
## H_box = scitools.BoxField.dolfin_function2BoxField(H_out, mesh_out, (M,N))
## H_values = H_box.values


nc.close()

try:
    move(tmp_filename, output_filename)
except:
    print "Error moving %s to %s. Exiting..." % (tmp_filename,
                                                 output_filename)
