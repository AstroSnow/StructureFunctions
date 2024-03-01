import numpy as np
import dedalus.public as d3
import random
import h5py
import logging
logger = logging.getLogger(__name__)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dedalus.extras.plot_tools import *

Lx, Ly = 1, 1
Nx, Ny = 64, 64
Reynolds = 5e3

dealias = 3/2
stop_sim_time = 20
timestepper = d3.RK443
max_timestep = 1e-5
dtype = np.float64

coords = d3.CartesianCoordinates('x', 'y')
dist = d3.Distributor(coords, dtype=dtype)
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(-Lx/2, Lx/2), dealias=dealias)
ybasis = d3.RealFourier(coords['y'], size=Ny, bounds=(-Ly/2, Ly/2), dealias=dealias)

A = dist.Field(name='A', bases=(xbasis,ybasis)) # in 3D
psi = dist.Field(name='psi', bases=(xbasis,ybasis)) # in 3D
c = dist.Field(name='c')

nu = 1 / Reynolds
eta = nu
x, y = dist.local_grids(xbasis, ybasis)
ex, ey = coords.unit_vector_fields(dist)

B = -d3.skew(d3.grad(A)) # dA/dx, dA/dy -> dA/dy, -dA/dx (d3. skew 90 degree positive rotation)
u = -d3.skew(d3.grad(psi))
w = -d3.lap(psi)
j = -d3.div(d3.skew(B)) # check if this is correct
problem = d3.IVP([A,psi,c], namespace=locals())
problem.add_equation("dt(w) - nu*lap(w) + c = - u@grad(w) + dot(B,grad(j))")
problem.add_equation("dt(A) - eta*lap(A)  = -dot(u,grad(A))")# u*grad(A) ")
problem.add_equation("integ(psi) = 0") # Pressure gauge



# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time



psi['g'] = np.sin(2*np.pi*x/Lx)
A['g'] = np.sin(2*np.pi*x/Lx)



# CFL
CFL = d3.CFL(solver, initial_dt=max_timestep, cadence=10, safety=0.2, threshold=0.1,
             max_change=1.5, min_change=0.5, max_dt=max_timestep)
CFL.add_velocity(u)

# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=10)
flow.add_property((u@ey)**2, name='w2')

# Main loop
try:
    logger.info('Starting main loop')
    while solver.proceed:
        timestep = CFL.compute_timestep()
        solver.step(timestep)
        if (solver.iteration-1) % 10 == 0:
            max_w = np.sqrt(flow.max('w2'))
            logger.info('Iteration=%i, Time=%e, dt=%e, max(w)=%f' %(solver.iteration, solver.sim_time, timestep, max_w))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()





