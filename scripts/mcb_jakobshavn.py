#!/usr/bin/env python

from dolfin import *


## TODO: Re-write to read in netCDF instead of Matlab mat file

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
