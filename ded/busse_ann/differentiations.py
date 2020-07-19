import numpy as np


def xdifferentiate(sslvr,fu,fieldnum=3):

    if fieldnum == 1:
        return xdiffOneField(sslvr,fu)
    elif fieldnum == 2:
        return xdiffTwoFields(sslvr,fu)
    elif fieldnum == 3:
        return xdiffThreeFields(sslvr,fu)
    else:
        raise NotImplementedError()

def xdiffThreeFields(solver,fu):

    domain = solver.domain

    psi   = solver.state['psi']
    theta = solver.state['theta']
    zeta = solver.state['zeta']

    psi.set_scales(1)
    theta.set_scales(1)
    zeta.set_scales(1)

    Nx, Ny = domain.local_grid_shape(scales=1)

    psi['g']   = fu[     0:Nx,0:Ny]
    theta['g'] = fu[  Nx:2*Nx,0:Ny]
    zeta['g']  = fu[2*Nx:3*Nx,0:Ny]


    cx, cy = domain.local_coeff_shape


    for k in range(cx):
        for i in range(cy):
            psi['c'][k , i]   = 1j *   psi['c'][k,i] *  domain.elements(0)[k,0]
            theta['c'][k ,i]  = 1j * theta['c'][k,i] *  domain.elements(0)[k,0] 
            zeta['c'][k , i]  = 1j * zeta['c'][k, i] *  domain.elements(0)[k,0]

    psix   = psi['g']
    thetax = theta['g']
    zetax  = zeta['g']

    result = np.vstack((psix, thetax, zetax))

    return result

def xdiffTwoFields(solver,fu):

    domain = solver.domain

    psi   = solver.state['psi']
    theta = solver.state['theta']
    zeta = solver.state['zeta']

    psi.set_scales(1)
    theta.set_scales(1)
    zeta.set_scales(1)

    Nx, Ny = domain.local_grid_shape(scales=1)

    psi['g']   = fu[     0:Nx,0:Ny]
    theta['g'] = fu[  Nx:2*Nx,0:Ny]
    zeta['g']  = fu[2*Nx:3*Nx,0:Ny]

    cx, cy = domain.local_coeff_shape

    for k in range(cx):
        for i in range(cy):
            psi['c'][k , i]   =   psi['c'][k,i] *  domain.elements(0)[k,0]
            theta['c'][k ,i]  = theta['c'][k,i] *  domain.elements(0)[k,0] 

    psix   = psi['g']
    thetax = theta['g']
    zet    = zeta['g']

    result = np.vstack((psix, thetax, zet))

    return result

def xdiffOneField(solver,fu):

    domain = solver.domain

    Nx, Ny = domain.local_grid_shape(scales=1)

    psi   = solver.state['psi']
    theta = solver.state['theta']
    zeta = solver.state['zeta']

    psi.set_scales(1)
    theta.set_scales(1)
    zeta.set_scales(1)

    psi['g']   = fu[0:Nx, 0:Ny]


    cx, cy = domain.local_coeff_shape


    for k in range(cx):
        for i in range(cy):
            psi['c'][k , i]   =   psi['c'][k,i] *  domain.elements(0)[k,0]

    psix   = psi['g']

    return psix
