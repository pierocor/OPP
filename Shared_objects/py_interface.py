#!/usr/bin/env python

from ctypes import *
testlib = CDLL('./py_shared.so')
natoms=10
class Data(Structure):
    _fields_ = [("natoms",c_int),("nfi", c_int),("nsteps", c_int),("epsilon",c_double),
               ("sigma",c_double),("dt",c_double),("mass",c_double),("box",c_double),("rcut",c_double),
               ("ekin",c_double),("epot",c_double),("temp",c_double),
               ("rx",POINTER(c_double)),("ry",POINTER(c_double)),("rz",POINTER(c_double)),
               ("vx",POINTER(c_double)),("vy",POINTER(c_double)),("vz",POINTER(c_double)),
               ("fx",POINTER(c_double)),("fy",POINTER(c_double)),("fz",POINTER(c_double))]

rx=(c_double * natoms)()
ry=(c_double * natoms)()
rz=(c_double * natoms)()
vx=(c_double * natoms)()
vy=(c_double * natoms)()
vz=(c_double * natoms)()
fx=(c_double * natoms)()
fy=(c_double * natoms)()
fz=(c_double * natoms)()
for i in range(natoms):
    rx[i] = i
    ry[i] = i
    rz[i] = i
    vx[i] = 0
    vy[i] = 0
    vz[i] = 0
    fx[i] = 0
    fy[i] = 0
    fz[i] = 0


system = Data( natoms = 10, nsteps = 10,epsilon = 0.2379, sigma = 3.405, dt = 5.0,\
 box = 17.158, rcut = 8.5, temp = 8.0, ekin=3.0, nfi=1, rx=rx, ry=ry, rz=rz, vx=vx, vy=vy, vz=vz, fx=fx, fy=fy, fz=fz)
testlib.output(byref(system))
