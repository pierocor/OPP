# -*- Makefile -*-
SHELL=/bin/sh
############################################
# derived makefile variables
OBJ_SERIAL=$(SRC:src/%.f90=Obj-serial/%.o)
############################################

default: mpi

serial:
	$(MAKE) $(MFLAGS) -C Obj-$@

mpi:
	$(MAKE) $(MFLAGS) -C Obj-$@

clean:
	$(MAKE) $(MFLAGS) -C Obj-serial clean
	$(MAKE) $(MFLAGS) -C Obj-mpi clean
	$(MAKE) $(MFLAGS) -C examples clean
	$(MAKE) $(MFLAGS) -C unit_test/force clean
	$(MAKE) $(MFLAGS) -C unit_test/kinetic clean
	$(MAKE) $(MFLAGS) -C unit_test/time_step clean
	$(MAKE) $(MFLAGS) -C unit_test/input clean


check: mpi default
	$(MAKE) $(MFLAGS) -C examples check

check_force:
	$(MAKE) $(MFLAGS) -C unit_test/force check

check_kinetic:
	$(MAKE) $(MFLAGS) -C unit_test/kinetic check

check_input:
	$(MAKE) $(MFLAGS) -C unit_test/input check

check_step:
	$(MAKE) $(MFLAGS) -C unit_test/time_step check
