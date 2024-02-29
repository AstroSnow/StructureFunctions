import numpy as np
import dedalus.public as d3
import random
import logging
logger = logging.getLogger(__name__)

# Parameters
Nx = 128
Ny = 128

Re = 1.0e4
Pr = 1.0e-1

dealias = 3/2
stop_sim_time = 1.0
timestepper = d3.RK443
max_timestep = 1e-2
dtype = np.float64

# Bases
coords = d3.CartesianCoordinates('x','y','z')
dist = d3.Distributor(coords, dtype=dtype)
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, 1), dealias=dealias)
ybasis = d3.RealFourier(coords['y'], size=Ny, bounds=(0, 1), dealias=dealias)

# Fields
p = dist.Field(name='p', bases=(xbasis,ybasis))
phi = dist.Field(name='phi', bases=(xbasis,ybasis))
u = dist.VectorField(coords, name='u', bases=(xbasis,ybasis))
A = dist.VectorField(coords, name='A', bases=(xbasis,ybasis))
tau_p = dist.Field(name='tau_p')
tau_phi = dist.Field(name='tau_phi')


# Substitutions
nu = 1/Re
eta = 1 / (Re*Pr)
x, y = dist.local_grids(xbasis, ybasis)
ex,ey,ez = coords.unit_vector_fields(dist)
dx = lambda A: d3.Differentiate(A, coords['x'])
dy = lambda A: d3.Differentiate(A, coords['y'])

#rho_h = dist.Field()
#rho_h['g'] = 3.0

#rho_l = dist.Field()
#rho_l['g'] = 1.0

a = dist.Field(bases = (xbasis,ybasis))
a2 = dist.Field(bases = (xbasis,ybasis))

B = d3.curl(A) # Is this not an equation?
problem = d3.IVP([u, A, phi, p, tau_p, tau_phi], namespace=locals())
problem.add_equation("dt(u) + grad(p) - nu*lap(u) = - u@grad(u) + B@grad(B)")
problem.add_equation("dt(A) - eta*lap(A) + grad(phi) = cross(u,B) ")
problem.add_equation("div(u) + tau_p = 0")
problem.add_equation("div(A) + tau_phi = 0") # Coulomb gauge
problem.add_equation("integ(p) = 0") # Pressure gauge
problem.add_equation("integ(phi) = 0") # Electric potential gauge

# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

# Initial conditions
# Background shear
#rho['g'] = 1 #- 2/2 * (np.tanh((z-((9*Lz)/20))/0.05) + 1) + 2/2 * (np.tanh(z/0.05) + 1)
#rho['g']=0.1 # + 1/2 * (np.tanh((z-0.5)/0.1) - np.tanh((z+0.5)/0.1))
u['g'][0] = -np.sin(2.0*np.pi*y)
u['g'][1] = np.sin(2.0*np.pi*x)
u['g'][2] = 0.0

b0=1.0
#b['g'][0] = -1.0/(4.0*np.pi)*np.sin(2.0*np.pi*y)
#b['g'][1] = 1.0/(4.0*np.pi)*np.sin(4.0*np.pi*x)

A['g'][0] = 0.0 #-1.0/(4.0*np.pi)*np.sin(2.0*np.pi*y)
A['g'][1] = 0.0 #1.0/(4.0*np.pi)*np.sin(4.0*np.pi*x)
A['g'][2] = 2.0*np.cos(2.0*np.pi*x)/(4.0*np.pi)**2+np.cos(4.0*np.pi*y)/(4.0*np.pi)**2


a['g'] = 0
a2['g'] = 0
for i in range(1,257):
	k = np.float128(i)
	amp = round(random.uniform(-1.0, 1.0), 5)
	phase = random.uniform(0, 2.0*np.pi)
	a['g'] = a['g'] + amp*np.sin((k*x) + phase)
	amp = round(random.uniform(-1, 1), 5)
	phase = random.uniform(0, 2.0*np.pi)
	a2['g'] = a2['g'] + amp*np.sin((k*x) + phase)
	
# Add small vertical velocity perturbations localized to the shear layers
u['g'][0] += a['g']*0.1
u['g'][1] += a2['g']*0.1
#u['g'][2] += 1.0


# Analysis
snapshots = solver.evaluator.add_file_handler('DataOZ', sim_dt=0.01, max_writes=1000)
#snapshots.add_task(rho, name='density')
#snapshots.add_task(d3.grad(rho), name='grad_density')
snapshots.add_task(p, name='pressure')
snapshots.add_task(u, name='u')
#snapshots.add_task(dx(u), name = 'dx_u')
#snapshots.add_task(dz(u), name = 'dz_u')
snapshots.add_task(B, name='B')
snapshots.add_task(A, name='A')
snapshots.add_task(phi, name='phi')
#snapshots.add_task(dx(b), name = 'dx_b')
#snapshots.add_task(dz(b), name = 'dz_b')
#snapshots.add_task(-d3.div(d3.skew(b)), name='j')
#snapshots.add_task(dx(-d3.div(d3.skew(b))), name='dx_j')
#snapshots.add_task(dz(-d3.div(d3.skew(b))), name='dz_j')
#snapshots.add_task(-d3.div(d3.skew(u)), name='vorticity')
#snapshots.add_task(4 * ((rho - rho_l)*(rho_h - rho))/((rho_h - rho_l)**2), name='theta') #Stone&Gardiner2007b,AstrophysicalJournal

# CFL
CFL = d3.CFL(solver, initial_dt=max_timestep*0.1, cadence=10, safety=0.2, threshold=0.1,
             max_change=1.5, min_change=0.5, max_dt=max_timestep)
CFL.add_velocity(u)

# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=10)
# Main loop
try:
    logger.info('Starting main loop')
    while solver.proceed:
        timestep = CFL.compute_timestep()
        solver.step(timestep)
        if (solver.iteration-1) % 10 == 0:
            logger.info('Iteration=%i, Time=%e, dt=%e' %(solver.iteration, solver.sim_time, timestep))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()

