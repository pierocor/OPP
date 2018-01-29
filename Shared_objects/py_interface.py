#!/usr/bin/env python
import ctypes
testlib = ctypes.cdll.LoadLibrary('./py_shared.so')
natoms=100
class Data(Structure):
    _filds_ = [("natoms",c_int),("nfi", c_int),("nsteps", c_int),("epsilon",c_float),
               ("sigma",c_float),("dt",c_float),("mass",c_float),("box",c_float),("rcut",c_float),
               ("ekin",c_float),("epot",c_float),("temp",c_float),
               ("rx",ctypes.POINTER(ctypes.c_float*natoms)),("ry",ctypes.POINTER(ctypes.c_float*natoms)),("rz",ctypes.POINTER(ctypes.c_float*natoms)),
               ("vx",ctypes.POINTER(ctypes.c_float*natoms)),("vy",ctypes.POINTER(ctypes.c_float*natoms)),("vz",ctypes.POINTER(ctypes.c_float*natoms)),
               ("fx",ctypes.POINTER(ctypes.c_float*natoms)),("fy",ctypes.POINTER(ctypes.c_float*natoms)),("fz",ctypes.POINTER(ctypes.c_float*natoms))]


system = Data( natoms = natoms, nsteps = 10,epsilon = 0.2379, sigma = 3.405, dt = 5.0, box = 17.158, rcut = 8.5)
testlib.verlet1(system)               