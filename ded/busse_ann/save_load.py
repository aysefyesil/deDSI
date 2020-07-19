import h5py
import numpy as np
from mpi4py import MPI
import os


def read_guess(ifilename, iteration):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank ==0:
        print("read_guess is called")

    it = iteration
    h5f = h5py.File(ifilename, 'r')

    psi   =  h5f['/tasks/psi'][it][:]
    theta =  h5f['/tasks/theta'][it][:]
    zeta  =  h5f['/tasks/zeta'][it][:]

    result = np.vstack((psi,theta,zeta))

    return result

def initiate_write_guess(folderpath="", filename=""):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank ==0:
        print("initiate_write_guess is called")
    if os.path.isdir(folderpath) !=True:
        os.mkdir(folderpath)

    ofilename = os.path.join(folderpath, filename+"_"+str(rank)+".h5")         

    with h5py.File(ofilename, "w") as h5f:
        tasks  = h5f.create_group('tasks')
        dset1 = tasks.create_dataset(name="psi",   dtype='float64', shape = (0, 0, 0), maxshape=(None, None, None))
        dset2 = tasks.create_dataset(name="theta", dtype='float64', shape = (0, 0, 0), maxshape=(None, None, None))
        dset3 = tasks.create_dataset(name="zeta",  dtype='float64', shape = (0, 0, 0), maxshape=(None, None, None))

    return ofilename


def append_guess(ifilename, guess):

    n, ny = np.shape(guess)

    nx = int(n/3)

    with h5py.File(ifilename, "a") as h5f:
        dset1 = h5f["tasks/psi"]
        dset2 = h5f["tasks/theta"]
        dset3 = h5f["tasks/zeta"]

        dset1.resize(dset1.shape[0] + 1,  axis=0)
        dset1.resize(nx, axis=1)
        dset1.resize(ny, axis=2)

        dset2.resize(dset2.shape[0] + 1,  axis=0)
        dset2.resize(nx, axis=1)
        dset2.resize(ny, axis=2)

        dset3.resize(dset3.shape[0] + 1,  axis=0)
        dset3.resize(nx, axis=1)
        dset3.resize(ny, axis=2)

        dset1[:] = guess[   0:nx,  :ny]
        dset2[:] = guess[  nx:2*nx,:ny]
        dset3[:] = guess[2*nx:3*nx,:ny]

    return 0

def write_guess(folderpath="", filename="", guess=0):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    n, ny = np.shape(guess)
    nx = int(n/3)
    
    comm.Barrier()
    if rank ==0:
        print("write_guess is called")
        if os.path.isdir(folderpath) !=True:
            os.mkdir(folderpath)
        import glob
        for filee in glob.glob(folderpath+"/"+filename+"_*[0-9]"):
            os.remove(filee) 

    comm.Barrier() 

    ofilename = os.path.join(folderpath, filename+"_"+str(rank)+".h5")

    with h5py.File(ofilename, "w") as h5f:
        tasks  = h5f.create_group('tasks')
        dset1 = tasks.create_dataset(name="psi",   dtype='float64', shape = (1, nx, ny))
        dset2 = tasks.create_dataset(name="theta", dtype='float64', shape = (1, nx, ny))
        dset3 = tasks.create_dataset(name="zeta",  dtype='float64', shape = (1, nx, ny))

        dset1[:] = guess[   0:nx,  :ny]
        dset2[:] = guess[  nx:2*nx,:ny]
        dset3[:] = guess[2*nx:3*nx,:ny]

    return 0



def save_guess(folderpath="./bestu", filename="best", guess=0):
    write_guess(folderpath, filename, guess)

      
    return 0
