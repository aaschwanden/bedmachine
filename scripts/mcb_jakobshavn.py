#!/usr/bin/env python

import numpy as np
from dolfin import *

try:
    from netCDF4 import Dataset as CDF
except:
    from netCDF3 import Dataset as CDF

try:
    import PyPISMTools.PyPISMTools as ppt
except:
    import PyPISMTools as ppt


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

output_order = ("x", "y")

data_file="jak_basin.nc"
nc_data = CDF(data_file, 'r')
xdim, ydim, zdim, tdim = get_dims(nc_vel)
rho = np.squeeze(permute(nc_data.variables["rho"], output_order=output_order))
H0 = np.squeeze(permute(nc_data.variables["thk"], output_order=output_order))
nc_data.close()

vel_file="jak_surf_vels.nc"
nc_vel = CDF(vel_file, 'r')
xdim, ydim, zdim, tdim = get_dims(nc_vel)
uvel = np.squeeze(permute(nc_vel.variables["us"], output_order=output_order))
vvel = np.squeeze(permute(nc_vel.variables["vs"], output_order=output_order))
nc_vel.close()

output_order = ("time", "x", "y")

bc_file="jak_input_v1.1.nc"
nc_bc = CDF(bc_file, 'r')
xdim, ydim, zdim, tdim = get_dims(nc_bc)
Hin = np.squeeze(permute(nc_bc.variables["thk"], output_order=output_order))
smb = np.squeeze(permute(nc_bc.variables["smb"], output_order=output_order))
nc_bc.close()

          
## from src.utilities import DataInput, DataOutput
## set_log_level(PROGRESS)

## # Import data
## dd = DataInput(1,"../data/Helheim/",\
##                 Mesh("../meshes/Helheim_mesh_500m.xml"),\
##                  "Helheim_SMB.mat","Helheim_Thickness_IDW.mat")
## H0  = dd.read("Helheim_Thickness_IDW.mat",flip=True)
## Hb  = dd.read("Helheim_Thickness_Bamber.mat",flip=True)
## S   = dd.read("Helheim_Surface_BPRC.mat",flip=True)
## u_o = dd.read("Helheim_vx.mat",flip=True)
## v_o = dd.read("Helheim_vy.mat",flip=True)
## adot= dd.read("Helheim_SMB.mat",flip=True)
## rho_d = dd.read("Helheim_rho.mat",flip=True,dg=True,bool_data=True)

# Prepare data for analysis
Hb.vector()[Hb.vector()<10.] = 10. # Lower limit on thickness

# Velocity norm
U = as_vector([u_o,v_o])  # unsmoothed
Unorm = project(sqrt(dot(U,U)) + 1e-10)

# Steepest descents velocities (if needed)
u_s = project(-S.dx(0) * Unorm)
v_s = project(-S.dx(1) * Unorm)

#u_threshold = 500000. # Norms of velocity greater than this are modeled 

#u_o.vector()[Unorm.vector()>u_threshold] = u_os.vector()[Unorm.vector()>u_threshold]
#v_o.vector()[Unorm.vector()>u_threshold] = v_os.vector()[Unorm.vector()>u_threshold]

# Velocity vector
#U = as_vector([u_o,v_o]) 
#Unorm = project(sqrt(dot(U,U)) + 1e-10)

# Ignore slow regions of ice.
utol = 5.0

def inside(x,on_boundary):
  return (Unorm(x[0],x[1]) < utol) or on_boundary 
   
dbc = DirichletBC(dd.func_space,Hb,inside)

# Solution and Trial function
H = Function(dd.func_space)
dH = TrialFunction(dd.func_space)

# Test Function
phi = TestFunction(dd.func_space)

# Misfit penalty.  Penalized differences between the calculated and observed thickness.

gamma = 50.0 # This is really large in this case because mass conservation is badly out.
# If this is lower, then the flight lines aren't obeyed near the terminus.
# This is something interesting to consider, mcb in cases where dS/dt > 100 m/a

# Regularization parameter (penalty on the gradient of the solution)
alpha = 2.

#Objective function
I = (div(U*H) - adot)**2*dx \
    + gamma*rho_d*0.5*(H-H0)**2*dx \
    + alpha*(H.dx(0)**2 + H.dx(1)**2)*dx

delta_I = derivative(I,H,phi)

J = derivative(delta_I,H,dH)

solve(delta_I==0,H,dbc,J=J)

# File ouput
do = DataOutput('../results/Helheim/')
data_out = {'mcb_bed':project(S-H),'bamber_bed':project(S-Hb),'flux_div':project(div(U*H)),'U':Unorm,'density':rho_d,\
        'u_o':u_o,'v_o':v_o,'adot':adot,'S':S,'H0':H0,'rho_d':rho_d,'thick':H,'delta_H':project((H-H0)*rho_d)}
do.write_dictionary_of_files(data_out)
