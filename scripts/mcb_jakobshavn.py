#!/usr/bin/env python

import numpy as np
from scipy.io import netcdf_file as CDF
import dolfin as dolf

def generate_expression_from_gridded_data(x, y, var):
    from scipy.interpolate import RectBivariateSpline
    interpolant = RectBivariateSpline(x, y, var)
    class newExpression(dolf.Expression):
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

# Note:
# We don't need to know the dimension order in the netCDF file
# because we get the dimension via get_dims, and then perumute
# to the desired order via permute.

# minimum ice thickness
thk_min = 10
output_order = ("x", "y")

data_file="jak_basin.nc"
nc_data = CDF(data_file, 'r')
xdim, ydim, zdim, tdim = get_dims(nc_data)
x = nc_data.variables[xdim][:]
y = nc_data.variables[ydim][:]
rho = np.squeeze(permute(nc_data.variables["rho"], output_order=output_order))
H0 = np.squeeze(permute(nc_data.variables["thk"], output_order=output_order))
M = len(x)
xmin = x[0]
xmax = x[-1]
N = len(y)
ymin = y[0]
ymax = y[-1]
nc_data.close()

# TODO: make fill value stuff more robust
vel_file="jak_surf_vels.nc"
nc_vel = CDF(vel_file, 'r')
xdim, ydim, zdim, tdim = get_dims(nc_vel)
uvel = np.squeeze(permute(nc_vel.variables["us"], output_order=output_order)).copy()
uvel_fill = nc_vel.variables["us"]._FillValue
uvel[uvel==uvel_fill] = 0.
vvel = np.squeeze(permute(nc_vel.variables["vs"], output_order=output_order)).copy()
vvel_fill = nc_vel.variables["vs"]._FillValue
vvel[uvel==vvel_fill] = 0.
nc_vel.close()

output_order = ("time", "x", "y")
bc_file="jak_input_v1.1.nc"
nc_bc = CDF(bc_file, 'r')
xdim, ydim, zdim, tdim = get_dims(nc_bc)
Hobs = np.squeeze(permute(nc_bc.variables["thk"], output_order=output_order)).copy()
Hobs[Hobs<thk_min] = thk_min
# No, we don't really want the upper surface fromm the SeaRISE data set.
# But it will do for testing puropses.
S = np.squeeze(permute(nc_bc.variables["usrf"], output_order=output_order))
smb = np.squeeze(permute(nc_bc.variables["smb"], output_order=output_order))
nc_bc.close()

dolf.set_log_level(dolf.PROGRESS)

scale_factor = 1
MM = int(np.ceil(M * scale_factor))
NN = int(np.ceil(N * scale_factor))
x_scaled = np.linspace(xmin, xmax, MM)
y_scaled = np.linspace(ymin, ymax, NN)
mesh = dolf.RectangleMesh(np.float(xmin), np.float(ymin),
                        np.float(xmax), np.float(ymax),
                        MM,
                        NN
                        )

dolf.parameters["form_compiler"]["optimize"] = True
func_space =  dolf.FunctionSpace(mesh,"CG", 1)
func_space_dg =  dolf.FunctionSpace(mesh,"DG", 1)

Hin = dolf.project(generate_expression_from_gridded_data(x, y, Hobs), func_space)
H0 = dolf.project(generate_expression_from_gridded_data(x, y, H0), func_space)
S = dolf.project(generate_expression_from_gridded_data(x, y, S), func_space)
adot = dolf.project(generate_expression_from_gridded_data(x, y, smb), func_space)
rho_d = dolf.project(generate_expression_from_gridded_data(x, y, rho), func_space_dg)
u_o = dolf.project(generate_expression_from_gridded_data(x, y, uvel), func_space)
v_o = dolf.project(generate_expression_from_gridded_data(x, y, vvel), func_space)

# Velocity norm
U = dolf.as_vector([u_o,v_o])  # unsmoothed
Unorm = dolf.project(dolf.sqrt(dolf.dot(U,U)) + 1e-10)

# Steepest descents velocities (if needed)
u_s = dolf.project(-S.dx(0) * Unorm)
v_s = dolf.project(-S.dx(1) * Unorm)

#u_threshold = 500000. # Norms of velocity greater than this are modeled 

#u_o.vector()[Unorm.vector()>u_threshold] = u_os.vector()[Unorm.vector()>u_threshold]
#v_o.vector()[Unorm.vector()>u_threshold] = v_os.vector()[Unorm.vector()>u_threshold]

# Velocity vector
#U = as_vector([u_o,v_o]) 
#Unorm = project(sqrt(dot(U,U)) + 1e-10)

# Ignore slow regions of ice.
utol = 5.0

def inside(x, on_boundary):
  return (Unorm(x[0],x[1]) < utol) or on_boundary 
   
dbc = dolf.DirichletBC(func_space, Hin, inside)

# Solution and Trial function
H = dolf.Function(func_space)
dH = dolf.TrialFunction(func_space)

# Test Function
phi = dolf.TestFunction(func_space)

# Misfit penalty.  Penalized differences between the calculated and observed thickness.

gamma = 50.0 # This is really large in this case because mass conservation is badly out.
# If this is lower, then the flight lines aren't obeyed near the terminus.
# This is something interesting to consider, mcb in cases where dS/dt > 100 m/a

# Regularization parameter (penalty on the gradient of the solution)
alpha = 2.

# Objective function
I = (dolf.div(U*H) - adot)**2*dolf.dx \
    + gamma*rho_d*0.5*(H-H0)**2*dolf.dx \
    + alpha*(H.dx(0)**2 + H.dx(1)**2)*dolf.dx

delta_I = dolf.derivative(I, H, phi)

J = dolf.derivative(delta_I, H, dH)

params = dolf.NonlinearVariationalSolver.default_parameters()
params['newton_solver']['maximum_iterations'] = 20

dolf.solve(delta_I==0, H, dbc, J=J, solver_parameters=params)
