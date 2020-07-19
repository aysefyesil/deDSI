"""
2D Busse Annulus
ref: Brummell & Hart (1993) fig 5a

"""
import numpy as np
from mpi4py import MPI

from dedalus import public as de
from dedalus.extras import flow_tools

import logging
from . import save_load as sl
from . import observables as ob

logger = logging.getLogger(__name__)


def set_mu(mu):
    global Ra
    Ra = mu

def get_mu():
    return Ra

# Create bases and domain
def build_solver(nx, ny, Lx, Ly, RA, BE, PR, CC):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    Nx = nx
    Ny = ny


    #for equilibrium:
    # C = 0, beta = 0, Pr = 1, Ra = 5000
    Ra = RA 
    Pr = PR
    beta = BE
    C = CC

    #for travelling wave
    #Ra = 7.542e5 Pr = 1. beta = 3E4 C = 0

    if rank == 0:
        print("dedalus: solver is built for: Nx "+str(nx)+" Ny: "+str(ny)+" Lx: "+str(Lx)+ " Ly: "+str(Ly)+" and Ra: "+str(Ra))
        print("dedalus: Pr "+str(Pr)+" beta: "+str(beta)+" C: "+str(C))

    x_basis = de.Fourier('x', Nx, interval=(0, Lx), dealias=3/2)
    y_basis = de.SinCos ('y', Ny, interval=(0, Ly), dealias=3/2)

    domain = de.Domain([x_basis, y_basis], grid_dtype=np.float64)


# 2D Boussinesq hydrodynamics
    problem = de.IVP(domain, variables=['psi','theta','zeta'], time='t')
    problem.meta['psi','zeta','theta']['y']['parity'] = -1 # sin basis
    problem.parameters['Pr'] = Pr
    problem.parameters['Ra'] = Ra
    problem.parameters['beta'] = beta
    problem.parameters['sbeta'] = np.sqrt(np.abs(beta))
    problem.parameters['C'] = C

# construct the 2D Jacobian
    problem.substitutions['J(A,B)'] = "dx(A) * dy(B) - dy(A) * dx(B)"

    problem.add_equation("dt(zeta) - beta*dx(psi) + Ra/Pr * dx(theta) + C * sbeta * zeta - dx(dx(zeta)) - dy(dy(zeta)) = -J(psi,zeta)", condition="ny != 0")
    problem.add_equation("dt(theta) + dx(psi) - (dx(dx(theta)) + dy(dy(theta)))/Pr = -J(psi,theta)", condition="ny != 0")
    problem.add_equation("dx(dx(psi)) + dy(dy(psi)) - zeta = 0", condition="ny != 0")
    problem.add_equation("zeta = 0", condition="ny ==0")
    problem.add_equation("theta = 0", condition="ny ==0")
    problem.add_equation("psi = 0", condition="ny ==0")

# Build solver
    solver = problem.build_solver(de.timesteppers.MCNAB2)
    #solver = problem.build_solver(de.timesteppers.SBDF1)
    logger.info('Solver built')
    return solver

def get_local_info(solver):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    nx, ny = solver.domain.local_grid_shape(scales=1)

    res = []

    res.append(int(rank))
    res.append(int(nx))
    res.append(int(ny))

    return res 

def distribute_guess(solver, ifilename, it, Amp=1):
    import h5py

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank ==0:
        print("distribute guess is called")

    slices = solver.domain.dist.grid_layout.slices(scales=1)

    h5f = h5py.File(ifilename, 'r')
    p =  h5f['/tasks/psi'][it][slices]
    t =  h5f['/tasks/theta'][it][slices]
    z =  h5f['/tasks/zeta'][it][slices]

    p = Amp * p
    t = Amp * t
    z = Amp * z
    result = np.vstack((p, t, z))

    return result

def distribute_guess2(solver, ifilename, it):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    uGuess = sl.read_guess(ifilename, it)


    NNx, Ny = np.shape(uGuess)

    Nx = int(NNx/3)

    uGuessp  = uGuess[:Nx, ]
    uGuesst  = uGuess[Nx:2*Nx, ]
    uGuessz  = uGuess[2*Nx:3*Nx, ]

    result = np.vstack((uGuessp, uGuesst, uGuessz))

    print(result)

    Nx, Ny = solver.domain.local_grid_shape(scales=1)

    psi    = solver.state['psi']
    theta  = solver.state['theta']
    zeta   = solver.state['zeta']

    psi.set_scales(1)
    theta.set_scales(1)
    zeta.set_scales(1)

    slices = solver.domain.dist.grid_layout.slices(scales=1)

    psi['g']   = uGuessp[slices]
    theta['g'] = uGuesst[slices]
    zeta['g']  = uGuessz[slices]

    result = np.vstack((psi['g'], theta['g'], zeta['g']))

    return result

def adjust_dt(T, dt):
    N = int(T/dt)
    last_dt = T - N * dt

    if last_dt == 0:
        last_dt = dt
        N = N - 1

    return N, last_dt

def integrate(solver, T , uGuess, dt_i):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    domain = solver.domain

    Nx, Ny = domain.local_grid_shape(scales=1)

    psi    = solver.state['psi']
    theta  = solver.state['theta']
    zeta   = solver.state['zeta']

    psi.set_scales(1)
    theta.set_scales(1)
    zeta.set_scales(1)

     
    psi['g']   = uGuess[:Nx, :Ny]
    theta['g'] = uGuess[Nx : 2*Nx, :Ny]
    zeta['g']  = uGuess[2*Nx : 3*Nx ,:Ny]
   

    solver.iteration = 0
    solver.timestepper._iteration = 0
    solver.sim_time  = 0
     
    logger.info("solver T: %.5e dt: %.8e "%(T,dt_i))
    
    N, last_dt = adjust_dt(T, dt_i)

    solver.stop_sim_time  = np.inf 
    solver.stop_wall_time = np.inf
    solver.stop_iteration = N + 1

    for handler in solver.evaluator.handlers:
        handler.last_wall_div = -1
        handler.last_sim_div  = -1
        handler.last_iter_div = -1

# Main loop

    dt = dt_i
    
    l2u = ob.l2norm(solver, uGuess)
    psi0 = psi['g'][0][0]

    logger.info('Iteration: %i, Time: %.17e, l2(u): %.17e, psi(0,0): %.17e dt: %.17e' %(solver.iteration, T, l2u, psi0, dt))
    
    try:
        logger.info('Starting loop')
        while solver.ok:
            if (solver.iteration <= N-1):
                dt = solver.step(dt)
            else:
                dt = solver.step(last_dt)

    except:
        logger.error('Exception raised, triggering end of main loop.')
        raise

    psi.set_scales(1)
    theta.set_scales(1)
    zeta.set_scales(1)

    result = np.vstack((psi['g'], theta['g'], zeta['g']))
    l2u = ob.l2norm(solver, result)
    psi0 = psi['g'][0][0]

    logger.info('Iteration: %i, Reached time: %.17e, expected time (calculated): %.17e target time: %.17e time difference: %.17e last_dt: %.17e dt: %.17e' %(solver.iteration, solver.sim_time, N * dt_i + last_dt , T, T-solver.sim_time, last_dt, dt_i))
    logger.info('l2(fu): %.17e, psi(0,0): %.17e ' %( l2u, psi0))

    return result

def normalization(solver, uGuess, norm_max, norm_l2, pnorm = 1.0, tnorm = 1.0, znorm = 1.0 ):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    Nx, Ny = solver.domain.local_grid_shape(scales=1)

    psi   = uGuess[   0:Nx  , 0:Ny]
    theta = uGuess[  Nx:2*Nx, 0:Ny]
    zeta  = uGuess[2*Nx:3*Nx, 0:Ny]

    #if (pnorm ==1.0 and tnorm == 1.0 and znorm == 1):
    if norm_max==True:
        pnorm1 = np.amax(np.abs(psi)) 
        tnorm1 = np.amax(np.abs(theta)) 
        znorm1 = np.amax(np.abs(zeta)) 
    
        pnorm = comm.allreduce(pnorm1, op=MPI.MAX) 
        tnorm = comm.allreduce(tnorm1, op=MPI.MAX) 
        znorm = comm.allreduce(znorm1, op=MPI.MAX) 
    
    elif norm_l2==True:
        pnorm = ob.l2norm(psi) 
        tnorm = ob.l2norm(theta) 
        znorm = ob.l2norm(zeta) 

    psi   /= pnorm 
    theta /= tnorm 
    zeta  /= znorm
    
    if rank == 0:
        if norm_max == True:
            print("dedalus normalization wrt maximum -> pnorm is "+str(pnorm)+" tnorm "+str(tnorm)+" znorm "+str(znorm)) 
        if norm_l2 == True:
            print("dedalus normalization wrt l2norm  -> pnorm is "+str(pnorm)+" tnorm "+str(tnorm)+" znorm "+str(znorm)) 

    result = np.vstack((psi, theta, zeta))

    return (result, pnorm, tnorm, znorm)

def unnormalization(solver, uGuess, pnorm = 1.0, tnorm = 1.0, znorm = 1.0 ):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    Nx, Ny = solver.domain.local_grid_shape(scales=1)

    psi   = uGuess[   0:Nx  , 0:Ny]
    theta = uGuess[  Nx:2*Nx, 0:Ny]
    zeta  = uGuess[2*Nx:3*Nx, 0:Ny]

    psi   *= pnorm 
    theta *= tnorm 
    zeta  *= znorm

    if rank == 0:
        print("dedalus unnormalization pnorm is "+str(pnorm)+" tnorm "+str(tnorm)+" znorm "+str(znorm)) 

    result = np.vstack((psi, theta, zeta))

    return result



