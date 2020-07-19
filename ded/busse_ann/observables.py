import numpy as np
from mpi4py import MPI
from . import differentiations as di

def energy(solver, guess):

    domain = solver.domain


    Nx, Ny = solver.domain.local_grid_shape(scales=1)

    start_int = solver.domain.bases[0].interval[0]
    end_int   = solver.domain.bases[0].interval[1]

    Lx = end_int - start_int

    start_y = domain.bases[1].interval[0]
    end_y   = domain.bases[1].interval[1]

    Ly = end_y - start_y

    psi   = solver.state['psi']
    psi.set_scales(1)

    psi['g'] = guess[0:Nx,0:Ny]

    psi['g'] = di.xdifferentiate(solver, psi['g'], 1)

    eng = 0

    p2_in = domain.new_field(name='p2_in')
    p2_in = psi**2
    p2_in = p2_in.evaluate()

    resp = 0
    resp = p2_in.integrate('x','y')
    resp = np.float(resp['g'][0][0]/(Lx*Ly))

    return resp 

def l2norm(solver, guess):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    domain = solver.domain

    start_x = domain.bases[0].interval[0]
    end_x   = domain.bases[0].interval[1]

    Lx = end_x - start_x

    start_y = domain.bases[1].interval[0]
    end_y   = domain.bases[1].interval[1]

    Ly = end_y - start_y

    Nx, Ny = solver.domain.local_grid_shape(scales=1)

    p = domain.new_field()
    p.meta['y']['parity'] = -1

    t = domain.new_field()
    t.meta['y']['parity'] = -1

    z = domain.new_field()
    z.meta['y']['parity'] = -1


    p.set_scales(1)
    t.set_scales(1)
    z.set_scales(1)

    p['g'] = guess[   0:Nx  , 0:Ny]
    t['g'] = guess[  Nx:2*Nx, 0:Ny]
    z['g'] = guess[2*Nx:3*Nx, 0:Ny]

    p2_in = p**2
    p2_in = p2_in.evaluate()
    p2_in = p2_in.integrate('x','y')
    resp = np.float(p2_in['g'][0][0])

    t2_in = t**2
    t2_in = t2_in.evaluate()
    t2_in = t2_in.integrate('x','y')
    rest = np.float(t2_in['g'][0][0])

    z2_in = z**2
    z2_in = z2_in.evaluate()
    z2_in = z2_in.integrate('x','y')
    resz = np.float(z2_in['g'][0][0])

    return np.sqrt((resp + rest + resz)/(Lx*Ly))


def extrimum_of_the_fields(solver, guess):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    Nx, Ny = solver.domain.local_grid_shape(scales=1)

    mx = np.zeros((3,1))
    mx_c = np.zeros((3,1))
    mn = np.zeros((3,1))
    mn_c = np.zeros((3,1))

    field1 = np.abs(guess[   0:Nx  , 0:Ny])
    field2 = np.abs(guess[  Nx:2*Nx, 0:Ny])
    field3 = np.abs(guess[2*Nx:3*Nx, 0:Ny])

    (mx[0,0], mn[0,0]) = (np.max(field1), np.min(field1))
    (mx[1,0], mn[1,0]) = (np.max(field2), np.min(field2))
    (mx[2,0], mn[2,0]) = (np.max(field3), np.min(field3))


    MPI.COMM_WORLD.Allreduce([mx, MPI.DOUBLE], [mx_c, MPI.DOUBLE], op=MPI.MAX)
    MPI.COMM_WORLD.Allreduce([mn, MPI.DOUBLE], [mn_c, MPI.DOUBLE], op=MPI.MIN)

    if rank == 0:
        print( " First field (psi) max: %.17e min: %.17e "%(mx_c[0,0],mn_c[0,0]))
        print( " Second field (theta) max: %.17e min: %.17e "%(mx_c[1,0],mn_c[1,0]))
        print( " Third field (zeta) max: %.17e min: %.17e "%(mx_c[2,0],mn_c[2,0]))

    return mx_c, mn_c


def l2norm_solver(solver, guess):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    domain = solver.domain

    start_x = domain.bases[0].interval[0]
    end_x   = domain.bases[0].interval[1]

    Lx = end_x - start_x

    start_y = domain.bases[1].interval[0]
    end_y   = domain.bases[1].interval[1]

    Ly = end_y - start_y

    Nx, Ny = solver.domain.local_grid_shape(scales=1)

    psi    = solver.state['psi']
    theta  = solver.state['theta']
    zeta   = solver.state['zeta']

    psi.set_scales(1)
    theta.set_scales(1)
    zeta.set_scales(1)

    psi['g']   = guess[   0:Nx  , 0:Ny]
    theta['g'] = guess[  Nx:2*Nx, 0:Ny]
    zeta['g']  = guess[2*Nx:3*Nx, 0:Ny]

    p2_in = psi**2
    p2_in = p2_in.evaluate()
    p2_in = p2_in.integrate('x','y')
    resp = np.float(p2_in['g'][0][0]/(Lx*Ly))
 
    t2_in = theta**2
    t2_in = t2_in.evaluate()
    t2_in = t2_in.integrate('x','y')
    rest = np.float(t2_in['g'][0][0]/(Lx*Ly))

    z2_in = zeta**2
    z2_in = z2_in.evaluate()
    z2_in = z2_in.integrate('x','y')
    resz = np.float(z2_in['g'][0][0]/(Lx*Ly))

    return np.sqrt((resp + rest + resz)/(Lx*Ly))


def l2norm_x1(Nx, Ny, guess):

    summ1 = np.zeros((1,1))
    gx = np.zeros((1,1))

    for i in range(3*Nx):
        for j in range(Ny):
            gx[0] += (guess[i][j])**2

    MPI.COMM_WORLD.Allreduce([gx, MPI.DOUBLE], [summ1, MPI.DOUBLE], op=MPI.SUM)

    return np.sqrt(summ1[0])[0]

def l2norm_x(Nx, Ny, guess):

    gx = np.array(0.0, 'd')
    summ1 = np.array(0.0, 'd')

    for i in range(3*Nx):
        for j in range(Ny):
            gx += guess[i][j]**2

    MPI.COMM_WORLD.Allreduce(gx, summ1, op=MPI.SUM)

    return np.sqrt(summ1)

def sorted_l2norm_x(guess):

    vxd = np.sort(guess, axis=None)    
    gx = np.array(0.0, 'd')
    summ1 = np.array(0.0, 'd')

    sz = np.size(vxd)
    print(sz)

    for i in range(sz):
        gx += vxd[i]**2
    MPI.COMM_WORLD.Allreduce(gx, summ1, op=MPI.SUM)

    return np.sqrt(summ1)

