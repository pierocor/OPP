# -*- Makefile -*-

default: check

clean:
	rm -f *.xyz *.dat

check: ../ljmd-serial.x
	../ljmd-serial.x < argon_108.inp
	head -10 argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> a.dat
	head -10 ../reference/argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> b.dat
	cmp a.dat b.dat || exit 1
	gprof ../ljmd-serial.x > ljmd-serial_argon_108.gprof   
	gprof ../ljmd-serial.x | ../gprof2dot-2017.9.19/gprof2dot.py | dot -T png -o callgraph_ljmd-serial-108.png 
	../ljmd-serial.x < argon_2916.inp
	head -10 argon_2916.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> a.dat
	head -10 ../reference/argon_2916.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> b.dat
	cmp a.dat b.dat || exit 1
	rm -f a.dat b.dat
	gprof ../ljmd-serial.x > ljmd-serial_argon_2916.gprof 
	gprof ../ljmd-serial.x | ../gprof2dot-2017.9.19/gprof2dot.py | dot -T png -o callgraph_ljmd-serial-2916.png 
