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
import scitools.BoxField

from argparse import ArgumentParser

def create_variable(var_name, f, long_name=None, standard_name=None,
                    units=None, fill_value=None, dimensions=None):
    var = nc.createVariable(var_name, np.double, dimensions)
    if long_name is not None:
        var.long_name = long_name
    if standard_name is not None:
        var.standard_name = standard_name
    if units is not None:
        var.units = units
    var.grid_mapping = "mapping"
    if fill_value is not None:
        var._FillValue = fill_value
    fb = scitools.BoxField.dolfin_function2BoxField(f, mesh, (MM,NN))
    var[0,:] = fb.values.T


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
parser.add_argument("-f","--scale_factor", dest="scale_factor",
                    type=float,
                    help='''Scale the element size Values<1 mean a coarser mesh. Default=1''',
                    default=1.)
parser.add_argument("--dhdt", dest="do_dhdt", action="store_true",
                  help="Include dhdt term.", default=False)
parser.add_argument("--bmelt", dest="do_bmelt", action="store_true",
                  help="Include basal melt term.", default=False)
parser.add_argument("-a","--alpha", dest="alpha", type=float,
                    help='''Regularization parameter''',
                    default=2.0)
parser.add_argument("-g","--gamma", dest="gamma", type=float,
                    help='''Misfit penalty''',
                    default=2.0)
parser.add_argument("-s","--velocity_scale", dest="velocity_scale", type=float,
                    help='''Scale the surface velocities''',
                    default=1.0)
parser.add_argument("--grid_spacing", dest="grid_spacing", type=int,
                    help='''Grid spacing in meters.''',
                    default=500)
parser.add_argument("-p","--project", dest="project_name",
                    help='''Name of the project which determines filenames together with grid size.''',
                    default='jakobshavn')

options = parser.parse_args()
scale_factor = options.scale_factor
do_dhdt = options.do_dhdt
do_bmelt = options.do_bmelt
parameters['allow_extrapolation'] = True
project_name = str(options.project_name)
grid_spacing = options.grid_spacing
# Misfit penalty.  Penalized differences between the calculated and observed thickness.
gamma = options.gamma
# Regularization parameter (penalty on the gradient of the solution)
alpha = options.alpha
# Velocity scale
velocity_scale = options.velocity_scale
# minimum ice thickness
thk_min = 10

output_order = ("x", "y")
filename = project_name + '_flightlines_' + str(grid_spacing) + 'm.nc'
nc = CDF(filename, 'r')
xdim, ydim, zdim, tdim = get_dims(nc)
x = nc.variables[xdim][:]
y = nc.variables[ydim][:]
proj4_str = nc.projection
input_dimensions = nc.variables["thk"].dimensions
rho = np.squeeze(permute(nc.variables["rho"], output_order=output_order))
Hfl = np.squeeze(permute(nc.variables["thk"], output_order=output_order))
M = len(x)
xmin = x[0]
xmax = x[-1]
N = len(y)
ymin = y[0]
ymax = y[-1]
nc.close()

output_order = ("x", "y")
filename = project_name + '_cresis_' + str(grid_spacing) + 'm.nc'
nc = CDF(filename, 'r')
xdim, ydim, zdim, tdim = get_dims(nc)
Hcresis = np.squeeze(permute(nc.variables["thk"], output_order=output_order))
Scresis = np.squeeze(permute(nc.variables["usurf"], output_order=output_order))
nc.close()

output_order = ("x", "y")
filename = project_name + '_surf_vels_' + str(grid_spacing) + 'm.nc'
nc = CDF(filename, 'r')
xdim, ydim, zdim, tdim = get_dims(nc)
uvel = np.squeeze(permute(nc.variables["us"], output_order=output_order))
try:
    uvel_fill = nc.variables["us"]._FillValue
    uvel[uvel.data==uvel_fill] = 0.
except:
    pass
#uvel[Scresis<thk_min] = 0.
vvel = np.squeeze(permute(nc.variables["vs"], output_order=output_order))
try:
    vvel_fill = nc.variables["vs"]._FillValue
    vvel[vvel.data==vvel_fill] = 0.
except:
    pass
#vvel[Scresis<thk_min] = 0.
nc.close()

output_order = ("x", "y")
filename = project_name + '_umt_' + str(grid_spacing) + 'm.nc'
nc = CDF(filename, 'r')
xdim, ydim, zdim, tdim = get_dims(nc)
Humt = np.squeeze(permute(nc.variables["thk"], output_order=output_order))
nc.close()

output_order = ("x", "y")
filename = project_name + '_bamber2013_' + str(grid_spacing) + 'm.nc'
nc = CDF(filename, 'r')
xdim, ydim, zdim, tdim = get_dims(nc)
Hbam13 = np.squeeze(permute(nc.variables["thk"], output_order=output_order))
nc.close()

output_order = ("time", "x", "y")
filename = project_name + '_searise_v1.1_' + str(grid_spacing) + 'm.nc'
nc = CDF(filename, 'r')
xdim, ydim, zdim, tdim = get_dims(nc)
Hsr = np.squeeze(permute(nc.variables["thk"], output_order=output_order))
smb = np.squeeze(permute(nc.variables["smb"], output_order=output_order))
nc.close()


output_order = ("x", "y")
filename = project_name + '_usurf_' + str(grid_spacing) + 'm.nc'
nc = CDF(filename, 'r')
xdim, ydim, zdim, tdim = get_dims(nc)
S = np.squeeze(permute(nc.variables["usurf"], output_order=output_order))
nc.close()

set_log_level(PROGRESS)

extend = 0
MM = int(np.ceil(M * scale_factor)) - 1
NN = int(np.ceil(N * scale_factor)) - 1
x_scaled = np.linspace(xmin - extend, xmax + extend, MM)
y_scaled = np.linspace(ymin - extend, ymax + extend, NN)
mesh = RectangleMesh(np.float(xmin) - extend, np.float(ymin - extend),
                        np.float(xmax + extend) ,np.float(ymax + extend),
                        MM, NN)

func_space =  FunctionSpace(mesh, "CG", 1)
func_space_dg =  FunctionSpace(mesh, "DG", 1)

H0_p = project(generate_expression_from_gridded_data(x, y, Hfl), func_space)
S_p = project(generate_expression_from_gridded_data(x, y, S), func_space)
smb_p = project(generate_expression_from_gridded_data(x, y, smb), func_space)
rho_p = project(generate_expression_from_gridded_data(x, y, rho, method="nearest"),
                func_space_dg)
Hcresis_p = project(generate_expression_from_gridded_data(x, y, Hcresis), func_space)
Hsr_p = project(generate_expression_from_gridded_data(x, y, Hsr), func_space)
Humt_p = project(generate_expression_from_gridded_data(x, y, Humt), func_space)
Hbam13_p = project(generate_expression_from_gridded_data(x, y, Hbam13), func_space)

u_o = project(generate_expression_from_gridded_data(x, y, uvel), func_space) * velocity_scale
v_o = project(generate_expression_from_gridded_data(x, y, vvel), func_space) * velocity_scale

if do_dhdt:
    output_order = ("x", "y")
    filename = project_name + '_dhdt_' + str(grid_spacing) + 'm.nc'
    nc = CDF(filename, 'r')
    xdim, ydim, zdim, tdim = get_dims(nc)
    dHdt = np.squeeze(permute(nc.variables["dhdt"], output_order=output_order))
    nc.close()
    dHdt_p = project(generate_expression_from_gridded_data(x, y, dHdt), func_space)
    dHdt_str = "dhdt"
else:
    dHdt_p = 0.
    dHdt_str = "nodhdt"
    
if do_bmelt:
    output_order = ("time", "x", "y")
    filename = project_name + '_bmelt_' + str(grid_spacing) + 'm.nc'
    nc = CDF(filename, 'r')
    xdim, ydim, zdim, tdim = get_dims(nc)
    bmelt = np.squeeze(permute(nc.variables["bmelt"], output_order=output_order))
    nc.close()
    bmelt_p = project(generate_expression_from_gridded_data(x, y, bmelt), func_space)
    bmelt_str = "bmelt"
else:
    bmelt_p = 0.
    bmelt_str = "nobmelt"

# Velocity norm
U = as_vector([u_o,v_o])  # unsmoothed
Unorm = project(sqrt(dot(U,U)) + 1e-10)

# Steepest descents velocities (if needed)
u_s = project(-S_p.dx(0) * Unorm)
v_s = project(-S_p.dx(1) * Unorm)

# Ignore slow regions of ice.
utol = 5.0

# Is this correct??
def inside(x, on_boundary):
  return (Unorm(x[0],x[1]) < utol)

def boundary(x, on_boundary):
  return (on_boundary and \
                       (x[0] < DOLFIN_EPS or x[1] < DOLFIN_EPS or \
                       (x[0] > 0.5 - DOLFIN_EPS and x[1] > 0.5 - DOLFIN_EPS)))

def bc_area(x, on_boundary):
  return (Unorm(x[0],x[1]) < utol) or (on_boundary and \
                       (x[0] < DOLFIN_EPS or x[1] < DOLFIN_EPS or \
                       (x[0] > 0.5 - DOLFIN_EPS and x[1] > 0.5 - DOLFIN_EPS)))

def bc_upstream(x, on_boundary):
    return on_boundary and (x[0] - 1 < DOLFIN_EPS)
            

## dbc_boundary = DirichletBC(func_space, Hcresis_p, boundary)
## dbc_inside = DirichletBC(func_space, 400., inside)
## dbc = [dbc_boundary, dbc_inside]

dbc_upstream = DirichletBC(func_space, Hcresis_p, bc_upstream)
dbc_area = DirichletBC(func_space, Hcresis_p, bc_area)
dbc = [dbc_area, dbc_upstream]

# Solution and Trial function
H = Function(func_space)
dH = TrialFunction(func_space)

# Test Function
phi = TestFunction(func_space)

# Objective function
I = 0.5*(div(U*H) - smb_p + bmelt_p + dHdt_p)**(2)*dx \
    + gamma*rho_p*0.5*(H-H0_p)**2*dx \
    + alpha*(H.dx(0)**2 + H.dx(1)**2)*dx

delta_I = derivative(I, H, phi)

J = derivative(delta_I, H, dH)

params = NonlinearVariationalSolver.default_parameters()
params['newton_solver']['relaxation_parameter'] = .6
params['newton_solver']['relative_tolerance'] = 1e-8
params['newton_solver']['absolute_tolerance'] = 1e-9
params['newton_solver']['maximum_iterations'] = 100

solve(delta_I==0, H, dbc, J=J, solver_parameters=params)

gamma_str = '_'.join(['gamma', str(gamma)])
alpha_str = '_'.join(['alpha', str(alpha)])
gs_str = str(grid_spacing) + 'm'
vel_str = '_'.join(['vscale', str(velocity_scale)])

import os
if not os.path.exists(project_name):
    os.makedirs(project_name)

out_dir = project_name + '/'
output_filename = '_'.join([project_name, gs_str, alpha_str,
                            gamma_str, vel_str, dHdt_str, bmelt_str,]) + '.nc'

mesh_out = RectangleMesh(np.float(xmin), np.float(ymin),
                         np.float(xmax), np.float(ymax),
                         M, N)


print "Creating output file..."
nc = CDF('/'.join([project_name, output_filename]), 'w')

nc.createDimension("y", size=y.shape[0])
nc.createDimension("x", size=x.shape[0])
nc.createDimension("time")

x_var = nc.createVariable("x", 'f', dimensions=("x",))
x_var.units = "m";
x_var.long_name = "easting"
x_var.standard_name = "projection_x_coordinate"
x_var[:] = x

y_var = nc.createVariable("y", 'f', dimensions=("y",))
y_var.units = "m";
y_var.long_name = "northing"
y_var.standard_name = "projection_y_coordinate"
y_var[:] = y

t_var = nc.createVariable("time", 'f', dimensions=("time",))
t_var.units = "years"
t_var.long_name = "time"
t_var.axis = "T"

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

dimensions = ("time", "y", "x")

create_variable("thk", project(H),
                long_name="land ice thickness from mass conservation",
                standard_name="land_ice_thickness",
                units="m", dimensions=dimensions)
create_variable("topg", project(S_p-H),
                long_name="bedrock surface elevation from mass conservation",
                standard_name="bedrock_altitude",
                units="m", dimensions=dimensions)
create_variable("U_mag", project(Unorm),
                long_name="magnitude of horizontal surface velocities",
                units="m year-1", dimensions=dimensions)
create_variable("divHU", project(div(H*U)),
                long_name="flux divergence using mass conservation",
                units="m year-1", dimensions=dimensions)
create_variable("divHU_cresis", project(div(Hcresis_p*U)), long_name="flux divergence cresis",
                units="m year-1", dimensions=dimensions)
create_variable("divHU_umt", project(div(Humt_p*U)), long_name="flux divergence UMT",
                units="m year-1", dimensions=dimensions)
create_variable("divHU_bamber2013", project(div(Hbam13_p*U)), long_name="flux divergence Bamber 2013",
                units="m year-1", dimensions=dimensions)
create_variable("divHU_searise", project(div(Hsr_p*U)), long_name="flux divergence SeaRISE",
                units="m year-1", dimensions=dimensions)
create_variable("divU", project(div(U)),
                long_name="divergence of velocity field",
                units="year-1", dimensions=dimensions)
create_variable("res_flux", project(div(U*H) - smb_p + bmelt_p + dHdt_p),
                long_name="residual flux",
                units="m2 year-1", dimensions=dimensions)
create_variable("cflux", project(sqrt((u_o*H)**2+(v_o*H)**2)),
                long_name="magnitude of vertically-averaged flux",
                units="m2 year-1", dimensions=dimensions)
create_variable("uflux", project(u_o*H),
                long_name="x-component of vertically-averaged flux",
                units="m2 year-1", dimensions=dimensions)
create_variable("vflux", project((v_o*H)),
                long_name="y-component of vertically-averaged flux",
                units="m2 year-1", dimensions=dimensions)

# Save the projection information:
nc.projection = proj4_str
nc.Conventions = "CF-1.5"

# writing global attributes
import time
va = vars(options)
script_command = ' '.join([time.ctime(), ':', __file__.split('/')[-1],
                           ' '.join([''.join(['--', str(l), ' ',str(va[l])]) for l in va])
])
nc.history = script_command

print "writing to %s ...\n" % '/'.join([project_name, output_filename])
print "run nc2cdo.py to add lat/lon variables" 
nc.close()


