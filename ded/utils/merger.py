import h5py
import os
import numpy as np
from mpi4py import MPI
import sys


_, folder_name = sys.argv

path_main   = os.getcwd()
path_folder = path_main+"/"+folder_name+"/"

filesindir = os.listdir(path_folder)

# sort the files
file_dict = {}
file_index = []

for fname in filesindir:
      file_index_value = int(fname.split('_')[1].split('.')[0])
      file_index.append(file_index_value)
      file_dict[file_index_value] = fname

lsorted = sorted(file_index)

sorted_files = []

for fidx in lsorted:
     sorted_files.append(file_dict[fidx])

# loop for getting into the folder and asking each file their size

g_ny = 0

for fname in sorted_files:
    with h5py.File(path_folder+fname,"r") as f:
        print("fname: ",f)
        g_it  = f['/tasks/psi'].shape[0]
        g_nx  = f['/tasks/psi'].shape[1]
        ny    = f['/tasks/psi'].shape[2]
        g_ny += ny 
        print("number of iterations: "+str(g_it)+" nx: "+str(g_nx)+" ny: "+str(ny)+" total ny: "+str(g_ny))


# creating the global output file with the size information gathered from the first step

fname_out = folder_name+"_all.h5"

with h5py.File(os.path.join(path_main, fname_out),"w") as f1:

    print("fname: ",f1)
    tasks = f1.create_group('tasks')
    dset1 = tasks.create_dataset(name="psi",   dtype='float64', shape = (g_it, g_nx, g_ny)) 
    dset2 = tasks.create_dataset(name="theta", dtype='float64', shape = (g_it, g_nx, g_ny))
    dset3 = tasks.create_dataset(name="zeta",  dtype='float64', shape = (g_it, g_nx, g_ny))

    ny_min = 0

    for fname in sorted_files:

        with h5py.File(path_folder+fname,"r") as f2:

            ny_max = ny_min + f2['/tasks/psi'].shape[2]

            print("fname: "+str(fname)+" ny_min: "+str(ny_min)+" ny_max: "+str(ny_max))

            dset1[:,:,ny_min:ny_max] = f2['/tasks/psi']
            dset2[:,:,ny_min:ny_max] = f2['/tasks/theta']
            dset3[:,:,ny_min:ny_max] = f2['/tasks/zeta']

            ny_min = ny_max

