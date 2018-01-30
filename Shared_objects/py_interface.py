#!/usr/bin/env python

from ctypes import *
testlib = CDLL('./py_shared.so')
natoms = 108
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

fh = open( "../examples/argon_108.inp" );

x = []
for line in fh.readlines():
    x.append( line.split()[0] )
fh.close()


system = Data()
system = Data( natoms = int(x[0]), mass = float(x[1]), nsteps = int(x[9]),epsilon = float(x[2]), sigma = float(x[3]),\
 dt = float(x[10]), box = float(x[5]), rcut = float(x[4]), temp = 8.0, ekin=3.0, nfi=1,  rx=rx, ry=ry, rz=rz, vx=vx, vy=vy, vz=vz, fx=fx, fy=fy, fz=fz)

nprint = int(x[11])

system.natoms = 3
system.nsteps = 10

print(system.natoms)

fh = open( "../examples/argon_108.rest" );

x1 = []
y1 = []
z1 = []
for line in fh.readlines():
	y = [value for value in line.split()]
	x1.append(float(y[0]))
	y1.append(float(y[1]))
	z1.append(float(y[2]))
fh.close()


for i in range(system.natoms):
	system.rx[i] = x1[i]
	system.ry[i] = y1[i]
	system.rz[i] = z1[i]
	system.vx[i] = 0.0 
	system.vy[i] = 0.0 
	system.vz[i] = 0.0
	system.fx[i] = 0.0
	system.fy[i] = 0.0 
	system.fz[i] = 0.0
	print (system.rx[i], system.ry[i], system.rz[i])

#system = Data( natoms = 10, nsteps = 10,epsilon = 0.2379, sigma = 3.405, dt = 5.0,\
# box = 17.158, rcut = 8.5, temp = 8.0, ekin=3.0, nfi=1, rx=rx, ry=ry, rz=rz, vx=vx, vy=vy, vz=vz, fx=fx, fy=fy, fz=fz)

#testlib.azzero(system.fx,system.natoms)
#testlib.azzero(system.fy,system.natoms)
#testlib.azzero(system.fz,system.natoms)


system.nfi=0

for i in range(system.natoms):
	system.fx[i] = 0.0
	system.fy[i] = 0.0 
	system.fz[i] = 0.0

testlib.force(byref(system))
testlib.ekin(byref(system))


#erg = open("output.erg","w")
#traj = open("output.traj","w")

for system.nfi in range(system.nsteps):
#	if (system.nfi+1) % nprint == 0:
#		print(system.nfi, system.temp, system.ekin, system.epot, system.ekin+system.epot)	
#		print(system.nfi, system.temp, system.ekin, system.epot, system.ekin+system.epot)
#		print(system.natoms, system.nfi,  system.ekin+system.epot)
#	print(2)
#	print(2)
	testlib.velverlet1(byref(system))
	for i in range(system.natoms):
		system.fx[i] = 0.0
		system.fy[i] = 0.0 
		system.fz[i] = 0.0	
	testlib.force(byref(system))
	testlib.velverlet2(byref(system))
	testlib.ekin(byref(system))

for i in range(system.natoms):
	print (system.rx[i], system.ry[i], system.rz[i])

print(3)

print("Simulation Done")	
#erg.close()
#traj.close()

print(4)
