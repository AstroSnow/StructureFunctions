import numpy as np
import dedalus.public as d3
import random
import logging
logger = logging.getLogger(__name__)

# Parameters
Lx, Ly = 10, 8
Nx = Lx*8	
Ny = Ly*8

Re = 1e2
Sc = 1	
Pr = 1

dealias = 3/2
stop_sim_time = 40
timestepper = d3.RK443
max_timestep = 1e-2
dtype = np.float64

# Bases
coords = d3.CartesianCoordinates('x','y', 'z')
dist = d3.Distributor(coords, dtype=dtype)
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(-Lx/2, Lx/2), dealias=dealias)
ybasis = d3.RealFourier(coords['y'], size=Ny, bounds=(-Ly/2, Ly/2), dealias=dealias)

# Fields
p = dist.Field(name='p', bases=(xbasis,ybasis))
rho = dist.Field(name='rho', bases=(xbasis,ybasis))
u = dist.VectorField(coords, name='u', bases=(xbasis,ybasis))
b = dist.VectorField(coords, name='b', bases=(xbasis,ybasis))
tau_u = dist.Field(name='tau_u')

# Substitutions
nu = 1 / Re
D = 1 / (Re*Sc)
eta = 1 / (Re*Pr)
x, y = dist.local_grids(xbasis, ybasis)
ex,ey, ez = coords.unit_vector_fields(dist)
dx = lambda A: d3.Differentiate(A, coords['x'])
dy = lambda A: d3.Differentiate(A, coords['y'])

#rho_h = dist.Field()
#rho_h['g'] = 3.0

#rho_l = dist.Field()
#rho_l['g'] = 1.0

#a = dist.Field(bases = (xbasis,zbasis))

# Problem

problem = d3.IVP([rho, u, tau_u, p, b], namespace=locals())
problem.add_equation("dt(rho) - D*lap(rho) = - u@grad(rho)")
problem.add_equation("dt(u) + grad(p) - nu*lap(u) = (b@grad(b))/rho - u@grad(u) - (p/rho)*grad(rho)")
problem.add_equation("div(u) + tau_u = 0")
problem.add_equation("dt(b) - eta*lap(b) - 100*grad(div(b)) = b@grad(u) - u@grad(b)")
problem.add_equation("integ(p) = 0")

# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

# Initial conditions
# Background shear
#rho['g'] = 1 #- 2/2 * (np.tanh((z-((9*Lz)/20))/0.05) + 1) + 2/2 * (np.tanh(z/0.05) + 1)
rho['g']=1 # + 1/2 * (np.tanh((z-0.5)/0.1) - np.tanh((z+0.5)/0.1))
u['g'][0] = 0.0 #1/2 + 1/2 * (np.tanh((z-0.5)/0.1) - np.tanh((z+0.5)/0.1))
u['g'][1] = 0.0
u['g'][2] = 0.0

b0=1.0
cw=1.0
yloc=2
xloc=5
wscale=4.0
bxpert=b0*x/2*np.exp(-(((x-xloc)**2+(y-yloc)**2)/(cw/wscale)))-b0*x/2*np.exp(-(((x+xloc)**2+(y+yloc)**2)/(cw/wscale)))
bypert=b0*y/2*np.exp(-(((x-xloc)**2+(y-yloc)**2)/(cw/wscale)))-b0*y/2*np.exp(-(((x+xloc)**2+(y+yloc)**2)/(cw/wscale)))
b['g'][0] = b0+b0*np.tanh((y-yloc)/cw)-b0*np.tanh((y+yloc)/cw)+bxpert
b['g'][1] = 0.0+bypert
b['g'][2] = b0/np.cosh((y-yloc)/cw)-b0/np.cosh((y+yloc)/cw)

"""a['g'] = 0
for i in range(1,257):
	k = np.float128(i)
	amp = round(random.uniform(-0.0125, 0.0125), 5)
	phase = random.uniform(0, np.pi)
	a['g'] = a['g'] + amp*np.sin((2*np.pi*k*x)/Lx + phase)
	
# Add small vertical velocity perturbations localized to the shear layers
u['g'][1] += a['g']*np.exp(-(z)**2/0.01) 
u['g'][2] = 1.0
"""

# Analysis
snapshots = solver.evaluator.add_file_handler('DataRecTest', sim_dt=0.5, max_writes=1000)
snapshots.add_task(rho, name='density')
snapshots.add_task(d3.grad(rho), name='grad_density')
snapshots.add_task(p, name='pressure')
snapshots.add_task(u, name='u')
#snapshots.add_task(dx(u), name = 'dx_u')
#snapshots.add_task(dz(u), name = 'dz_u')
snapshots.add_task(b, name='b')
#snapshots.add_task(dx(b), name = 'dx_b')
#snapshots.add_task(dz(b), name = 'dz_b')
#snapshots.add_task(-d3.div(d3.skew(b)), name='j')
#snapshots.add_task(dx(-d3.div(d3.skew(b))), name='dx_j')
#snapshots.add_task(dz(-d3.div(d3.skew(b))), name='dz_j')
#snapshots.add_task(-d3.div(d3.skew(u)), name='vorticity')
#snapshots.add_task(4 * ((rho - rho_l)*(rho_h - rho))/((rho_h - rho_l)**2), name='theta') #Stone&Gardiner2007b,AstrophysicalJournal

# CFL
CFL = d3.CFL(solver, initial_dt=max_timestep, cadence=10, safety=0.2, threshold=0.1,
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

