# -*- Makefile -*-
SHELL=/bin/sh
############################################
# derived makefile variables
OBJ_SERIAL=$(SRC:src/%.f90=Obj-serial/%.o)
############################################

default: serialx

serial:
	$(MAKE) $(MFLAGS) -C Obj-$@

mpi:
	$(MAKE) $(MFLAGS) -C Obj-parallel mpi

omp:
	$(MAKE) $(MFLAGS) -C Obj-parallel omp

clean:
	$(MAKE) $(MFLAGS) -C Obj-serial clean
	$(MAKE) $(MFLAGS) -C Obj-parallel clean
	$(MAKE) $(MFLAGS) -C examples clean
	$(MAKE) $(MFLAGS) -C unit_test/force clean
	$(MAKE) $(MFLAGS) -C unit_test/kinetic clean
	$(MAKE) $(MFLAGS) -C unit_test/time_step clean
	$(MAKE) $(MFLAGS) -C unit_test/input clean
	$(MAKE) $(MFLAGS) -C Shared_objects clean

py_wrapper:
	$(MAKE) $(MFLAGS) -C Shared_objects shared

py_check:
	$(MAKE) $(MFLAGS) -C Shared_objects check

check: serial
	$(MAKE) $(MFLAGS) -C examples check

check_mpi: mpi
	$(MAKE) $(MFLAGS) -C examples check_mpi

check_omp: omp
	$(MAKE) $(MFLAGS) -C examples check_omp

check_force:
	$(MAKE) $(MFLAGS) -C unit_test/force check

check_kinetic:
	$(MAKE) $(MFLAGS) -C unit_test/kinetic check

check_input:
	$(MAKE) $(MFLAGS) -C unit_test/input check

check_step:
	$(MAKE) $(MFLAGS) -C unit_test/time_step check
