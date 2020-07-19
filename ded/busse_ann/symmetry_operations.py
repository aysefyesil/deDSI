import numpy as np
#import mpi4py as MPI


from mpi4py import MPI

# sigma = [ s, sx, sy, ax ]

#def translation(domain, sigma, uguess):
#
#    ax = sigma[3]
#
#    start_int = domain.bases[0].interval[0]
#    end_int = domain.bases[0].interval[1]
#
#    L = end_int - start_int
#
#    v = domain.new_field()
#    v.set_scales(1)
#    v['g'] = uguess
#    print("layout:",v._layout)
#    print(v.atoms())
#    print("layout:",v._layout)
##    v.towards_grid_space
#    v['g'] 
#    print("layout:",v._layout)
#    print(type(v))
##    coeff = v.layout(-1)
#    print(v['g'])
#    v['c']
#
#    for k_x in range(len(domain.bases[0].elements)):
#        v['c'][k_x][:] = v['c'][k_x][:] * np.exp(1j * domain.bases[0].elements[k_x] * ax * L)
#
#    return v['g']
# symmetry transformations, for searching:
# the sigma notation is taken from channelflow see Channelflow wiki
# Here vorticity V = -dV/dy * e_x + dV/dx * e_y = u * e_x + v * e_y
# sigma [u,v] (x,y) = s( sx*u , sy*v)(sx * x + ax * Lx , sy * y)
# sigma [u,v] (x,y) = s( sx*u , sy*v)(sx * x + ax * Lx , sy * y)
# particularly for this problem possible symmetry transformations of the system are 
# J. Fluid Mech. (2006), vol. 567, pp. 117-140 section 4.3
# Busse annulus equation has two symmetries: 
# 1) [V,T](x,y) = [-V,T](-x,y) : V changes sign and x is reflected around x = Lx/2 and also scalar parameter beta -> -beta
# 2) [V,T](x,y) = -[V,T](x,-y) : both V and T change sign and y is reflected around y = Ly/2 
# And sigma[V,T](x,y)= s[ sx * V, T](sx*x + ax * Lx , sy * y) 
# sigma = [ s, sx, sy, ay]
# sigma1 = [1, -1, 1, ay != 0] provided that beta -> -beta
# sigma2 = [-1, 1, 1, -1, ay = 0] 
# Moreover since there is no y symmetry of the system in the first case,
# equations may support traveling wave solutions ini periodic x-direction 

# sigma = [ s, sx, sy, ax]

def translation(solver, ax, uguess, fieldnum=3):
    if fieldnum == 1:
       return translateSingleField(solver, ax, uguess)
    elif fieldnum == 2:
       return  translateDoubleField(solver, ax, uguess)
    elif fieldnum == 3:
       return  translateTripleField(solver, ax, uguess)
    else:
        raise NotImplementedError()

def translateSingleField(solver, sigma, uguess):
# sigma = [ s, sx, sy, ax ]
    ax = sigma[3]

    domain = solver.domain

    start_int = domain.bases[0].interval[0]
    end_int   = domain.bases[0].interval[1]

    L = end_int - start_int

    psi   = solver.state['psi']
    Nx, Ny = domain.local_grid_shape(scales=1)

    psi['g'] = uguess

    # add something which only allows translation operation in Fourier

    for k in range(len(domain.bases[0].elements)):
        psi['c'][k][:] = psi['c'][k][:] * np.exp(1j * domain.bases[0].elements[k] * ax * L)

    return psi['g']

def translateDoubleField(sslver, ax, uguess):

    solver = sslver
    domain = solver.domain

    Nx = len(domain.grid(0))

    start_int = domain.bases[0].interval[0]
    end_int   = domain.bases[0].interval[1]

    L = end_int - start_int

    psi   = solver.state['psi']
    theta = solver.state['theta']

    Nx, Ny = domain.local_grid_shape(scales=1)

    psi['g']   = uguess[0:Nx,:Ny]
    theta['g'] = uguess[Nx:2*Nx][:Ny]

    with open("psi.txt", 'w') as out:
        for i in range(Nx):
            for j in range(Ny):
                out.write("{} {} {} \n".format(i,j,psi['g'][i][j]))
            out.write("\n")
    out.close()

    with open("theta.txt", 'w') as fout:
        for i in range(Nx):
            for j in range(Ny):
                fout.write("{} {} {} \n".format(i,j,theta['g'][i][j]))
            fout.write("\n")
    fout.close()

    for k in range(len(domain.bases[0].elements)):
        psi['c'][k][:]   = psi['c'][k][:]   * np.exp(1j * domain.bases[0].elements[k] * ax * L)
        theta['c'][k][:] = theta['c'][k][:] * np.exp(1j * domain.bases[0].elements[k] * ax * L)

    psix   = psi['g']
    thetax = theta['g']

    result = np.vstack((psix, thetax))

    return result

def translateTripleField(sslver, ax, uguess):

    solver = sslver
    domain = solver.domain
    
    start_int = domain.bases[0].interval[0]
    end_int   = domain.bases[0].interval[1]

    L = end_int - start_int

    psi   = solver.state['psi']
    theta = solver.state['theta']
    zeta  = solver.state['zeta']

    psi.set_scales(1)
    theta.set_scales(1)
    zeta.set_scales(1)

    Nx, Ny = domain.local_grid_shape(scales=1)


    psi['g']   = uguess[   0:Nx,   :Ny]
    theta['g'] = uguess[  Nx:2*Nx, :Ny]
    zeta['g']  = uguess[2*Nx:3*Nx, :Ny]
    
    
    cx, cy = domain.local_coeff_shape
    
    for k in range(cx):
        for i in range(cy):
            psi['c'][k , i]   =   psi['c'][k,i] * np.exp(1j * domain.elements(0)[k,0] * ax * L)
            theta['c'][k ,i]  = theta['c'][k,i] * np.exp(1j * domain.elements(0)[k,0] * ax * L)
            zeta['c'][k , i]  = zeta['c'][k, i] * np.exp(1j * domain.elements(0)[k,0] * ax * L)

    result = np.vstack((psi['g'], theta['g'], zeta['g']))

    return result

def reflection(sslver, uguess, sp=1, st=1, sx=1, sy=1, case=None):
     if case==1:
         sp =-1
         st=1
         sx=-1
         sy=1

     elif case==2:
         sp=-1
         st=-1
         sx=1
         sy=-1

     print("sp is: ", sp," st is: ", st, " sx is: ", sx, " sy is: ", sy)

     solver = sslver
     Nx = len(solver.domain.grid(0))

     psi   = np.copy(uguess[:Nx][:])
     theta = np.copy(uguess[Nx:][:])
         

     if sx == -1: #flipping around (Lx/2) on the horizontal axis x
         psitemp = np.flipud(psi)
         psi = np.copy(psitemp)
         thetatemp = np.flipud(theta)
         theta = np.copy(thetatemp)
 
     if sy == -1: #flipping around (Ly/2) on the vertical axis y
         psitemp = np.fliplr(psi)
         psi = np.copy(psitemp)
         thetatemp = np.fliplr(theta)
         theta = np.copy(thetatemp)
     
     if sp == -1:
         psi = -psi

     if st == -1:
         theta = -theta
        
     result = np.vstack((psi, theta)) 

     return result 




