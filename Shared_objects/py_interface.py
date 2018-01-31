#!/usr/bin/env python

from ctypes import *
import time

testlib = CDLL('./py_shared.so')
natoms = 108
class Data(Structure):
    _fields_ = [("natoms",c_int),("nfi", c_int),("nsteps", c_int),("dt",c_double),("mass",c_double),("epsilon",c_double),
               ("sigma",c_double),("box",c_double),("rcut",c_double),
               ("ekin",c_double),("epot",c_double),("temp",c_double),
               ("rx",POINTER(c_double)),("ry",POINTER(c_double)),("rz",POINTER(c_double)),
               ("vx",POINTER(c_double)),("vy",POINTER(c_double)),("vz",POINTER(c_double)),
               ("fx",POINTER(c_double)),("fy",POINTER(c_double)),("fz",POINTER(c_double))]



def output( system, erg, traj):
    print(system.nfi, system.temp, system.ekin, system.epot, system.ekin+system.epot)
    erg.write("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n" % (system.nfi, system.temp, system.ekin, system.epot, system.ekin+system.epot))
    traj.write("%d\n nfi=%d etot=%20.8f\n" % (system.natoms, system.nfi, system.ekin+system.epot))
    for i in range(system.natoms):
        traj.write("Ar  %20.8f %20.8f %20.8f\n" % (system.rx[i], system.ry[i], system.rz[i]))

rx=(c_double * natoms)()
ry=(c_double * natoms)()
rz=(c_double * natoms)()
vx=(c_double * natoms)()
vy=(c_double * natoms)()
vz=(c_double * natoms)()
fx=(c_double * natoms)()
fy=(c_double * natoms)()
fz=(c_double * natoms)()

fh = open( "../examples/argon_108.inp" );

x = []
for line in fh.readlines():
    x.append( line.split()[0] )
fh.close()

system = Data()
system = Data( natoms = int(x[0]), mass = float(x[1]), nsteps = int(x[9]),epsilon = float(x[2]), sigma = float(x[3]),\
 dt = float(x[10]), box = float(x[5]), rcut = float(x[4]), rx=rx, ry=ry, rz=rz, vx=vx, vy=vy, vz=vz, fx=fx, fy=fy, fz=fz)

nprint = int(x[11])

fh = open( "../examples/argon_108.rest" );


for i in range(system.natoms):
    line=fh.readline()
    y = line.split()
    system.rx[i]=float(y[0])
    system.ry[i]=float(y[1])
    system.rz[i]=float(y[2])

for i in range(system.natoms):
    line=fh.readline()
    y = line.split()
    system.vx[i]=float(y[0])
    system.vy[i]=float(y[1])
    system.vz[i]=float(y[2])
fh.close()

erg = open(x[8], "w")
traj = open(x[7], "a")

for i in range(system.natoms):
	system.fx[i] = 0.0
	system.fy[i] = 0.0
	system.fz[i] = 0.0

testlib.force(byref(system))
testlib.ekin(byref(system))


print("Starting simulation with ",system.natoms," atoms for ",system.nsteps," steps." )
print("     NFI            TEMP            EKIN                 EPOT              ETOT")

t0 = time.time()
for system.nfi in range(1,system.nsteps+1):
    if (system.nfi % nprint == 0):
        output(system, erg, traj)
    for i in range(system.natoms):
        system.fx[i] = 0.0
        system.fy[i] = 0.0
        system.fz[i] = 0.0
    testlib.velverlet1(byref(system))
    testlib.force(byref(system))
    testlib.velverlet2(byref(system))
    testlib.ekin(byref(system))
t1 = time.time()


for i in range(system.natoms):
    print(system.rx[i], system.ry[i], system.rz[i])
    print(system.vx[i], system.vy[i], system.vz[i])
    print(system.fx[i], system.fy[i], system.fz[i])
print('Time taken:%20.8f sec \n' %(t1-t0))

print(system.natoms, system.mass, system.nsteps ,system.epsilon, system.sigma,\
 system.dt, system.box , system.rcut , system.temp, system.ekin, system.nfi)

erg.close()
traj.close()
print("Simulation Done")
