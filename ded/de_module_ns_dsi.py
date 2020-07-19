'''
Author: Ayse Yesil
License: GPL v2

This is an interface between channelflow's  nsolver and dedalus.

Below is the required list of functions
that is needed for nsolver interface

'''
import numpy as np
from mpi4py import MPI
import os
import sys
sys.path.append("./")

from busse_ann import time_integration as ti
from busse_ann import observables as ob
from busse_ann import symmetry_operations as so
from busse_ann import differentiations as di
from busse_ann import save_load as sl
## setting the continuation parameter mu:

def set_mu(mu):
    ti.set_mu(mu)
#
def get_mu():
    return ti.get_mu()
#
def build_solver(nx, ny, Lx, Ly, RA, BE, PR, CC):
    return ti.build_solver(nx, ny, Lx, Ly, RA, BE, PR, CC )

#returns information related to the paralelization
# such as local size and rank
def get_local_info(sslver):
    return ti.get_local_info(sslver)

def distribute_guess(solver, ifilename,it, amp):
    return ti.distribute_guess(solver, ifilename,it, amp)

def distribute_guess2(solver, ifilename,it):
    return ti.distribute_guess2(solver, ifilename,it)

def adjust_dt(T, dt):
    return ti.adjust_dt(T, dt)

#  u -> f(u) where f is the time integration of u (uGuess) upto T = n_iter * dt_i
def integrate(sslver, T, uGuess, dt_i):
    return ti.integrate(sslver, T, uGuess, dt_i)
#
def energy(sslver, uGuess):
    return ob.energy(sslver,uGuess)
#
def l2norm(solver, uGuess):
    return ob.l2norm(solver, uGuess)

def l2norm_x (Nx, Ny, guess):
    return ob.l2norm_x(Nx, Ny, guess)

## The following functions translation, reflection
## are only neccessary for non-equilibrium solutions (periodic orbits, travelling waves etc)

def translation(solver, ax, uGuess, fieldnum=3):
    # ax != 0
    return so.translation(solver, ax, uGuess, fieldnum=3)
#
def reflection(sslver, uguess, sp=1, st=1, sx=1, sy=1, case=None):
    #same ideas in the channelflow symmetry is at play here.
    # either s = -1 and one or none of sx, sy equal to -1
    # or s = 1 and one or two of sx, sy equal to -1
    return so.reflection(sslver, uguess, sp, st, sx, sy, case)
#
def xdifferentiate(sslvr, uGuess, fielnum):
    return di.xdifferentiate(sslvr, uGuess, fielnum)

def read_guess(ifilename, iteration):
    return sl.read_guess(ifilename, iteration)

def initiate_write_guess(ofoldername, ofilename):
    return sl.initiate_write_guess(ofoldername, ofilename)

def save_guess(ofoldername, ofilename, guess):
    return sl.save_guess(ofoldername, ofilename, guess)

def append_guess(ifilename, guess):
    return sl.append_guess(ifilename, guess)

def normalization(solver, uGuess, unorm_max, unorm_l2, pnorm, tnorm, znorm ):
    return ti.normalization(solver, uGuess, unorm_max, unorm_l2, pnorm, tnorm, znorm )

def unnormalization(solver, uGuess, pnorm, tnorm, znorm ):
    return ti.unnormalization(solver, uGuess, pnorm, tnorm, znorm )

def extrimum_of_the_fields(solver, guess):
    return ob.extrimum_of_the_fields(solver, guess)

def sorted_l2norm_x(guess):
    return ob.sorted_l2norm_x(guess)

def adjust_for_T(T, dt):
    n = round(T/dt)
    dt = T/n
    return dt
