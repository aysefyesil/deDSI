import numpy as np
from mpi4py import MPI


def display(a):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nx, ny = np.shape(a)
    if rank == 0:
        for i in range(nx):
            for j in range(ny):
                print("%d %d %.16e "%(i,j, a[i][j]))
    return

