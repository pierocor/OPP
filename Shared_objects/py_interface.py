#!/usr/bin/env python

from ctypes import *
testlib = CDLL('./py_shared.so')
class Data(Structure):
    _fields_ = [("natoms",c_int),("nfi", c_int),("nsteps", c_int),("epsilon",c_double),
               ("sigma",c_double),("dt",c_double),("mass",c_double),("box",c_double),("rcut",c_double),
               ("ekin",c_double),("epot",c_double),("temp",c_double),
               ("rx",POINTER(c_double)),("ry",POINTER(c_double)),("rz",POINTER(c_double)),
               ("vx",POINTER(c_double)),("vy",POINTER(c_double)),("vz",POINTER(c_double)),
               ("fx",POINTER(c_double)),("fy",POINTER(c_double)),("fz",POINTER(c_double))]


fh = open( "../examples/argon_108.inp" );

x = []
for line in fh.readlines():
    x.append( line.split()[0] )
fh.close()
print (x)

system = Data()
system = Data( natoms = int(x[0]), mass = float(x[1]), nsteps = int(x[9]),epsilon = float(x[2]), sigma = float(x[3]), dt = float(x[10]),\
 box = float(x[5]), rcut = float(x[4]), temp = 8.0, ekin=3.0, nfi=1, nprint = int(x[11]))

rx=(c_double * system.natoms)()
ry=(c_double * system.natoms)()
rz=(c_double * system.natoms)()
vx=(c_double * system.natoms)()
vy=(c_double * system.natoms)()
vz=(c_double * system.natoms)()
fx=(c_double * system.natoms)()
fy=(c_double * system.natoms)()
fz=(c_double * system.natoms)()

print(system.natoms)

fh = open( "../examples/argon_108.rest" );

x1 = []
y1 = []
z1 = []
for line in fh.readlines():
	y = [value for value in line.split()]
	x1.append(float(y[0]))
	y1.append(float(y[0]))
	z1.append(float(y[0]))
fh.close()



for i in range(system.natoms):
	rx[i] = x1[i]
	ry[i] = y1[i]
	rz[i] = z1[i]
#	print rx[i], ry[i], rz[i]
    
system = Data( natoms = 10, nsteps = 10,epsilon = 0.2379, sigma = 3.405, dt = 5.0,\
 box = 17.158, rcut = 8.5, temp = 8.0, ekin=3.0, nfi=1, rx=rx, ry=ry, rz=rz, vx=vx, vy=vy, vz=vz, fx=fx, fy=fy, fz=fz)
#testlib.output(byref(system))
testlib.velverlet1(byref(system))
testlib.velverlet2(byref(system))
testlib.force(byref(system))
testlib.output(byref(system))
