/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <data.h>
#include <velocity_verlet.h>

#ifdef MPI
#include <mpi.h>
#define TAG 100
#endif

/* a few physical constants */
extern const double kboltz;
extern const double mvsq2e;

/* main */
int main(int argc, char **argv)
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    double t_start, t_stop;
    mdsys_t sys;

#ifdef MPI
MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &sys.size);
MPI_Comm_rank(MPI_COMM_WORLD, &sys.rank);


/* INPUT PARAMETER */
if (  sys.rank == 0 ){
#endif

    /* read input file */
    if(get_a_line(stdin,line)) return 1;
    sys.natoms=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.mass=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.epsilon=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.sigma=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.rcut=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.box=atof(line);
    if(get_a_line(stdin,restfile)) return 1;
    if(get_a_line(stdin,trajfile)) return 1;
    if(get_a_line(stdin,ergfile)) return 1;
    if(get_a_line(stdin,line)) return 1;
    sys.nsteps=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.dt=atof(line);
    if(get_a_line(stdin,line)) return 1;
    nprint=atoi(line);

#ifdef MPI
}
/* BROADCAST PARAMETER */
MPI_Bcast( &sys.natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast( &sys.nfi, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast( &sys.nsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast( &sys.dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast( &sys.mass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast( &sys.epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast( &sys.sigma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast( &sys.box, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast( &sys.rcut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


/* COMPUTING RANGES */
// the standard range of a section is N/ sys.size
// if N is not a multiple of  sys.size the first processes will have an extra atom
sys.my_range = sys.natoms /  sys.size;
if( sys.rank < (sys.natoms %  sys.size) )
  sys.my_range++;
sys.my_start = sys.rank * (sys.natoms / sys.size);
if ( sys.rank < (sys.natoms % sys.size) )
  sys.my_start += sys.rank;
else
  sys.my_start += (sys.natoms % sys.size);
sys.my_end = sys.my_start + sys.my_range;
#endif

    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
#ifdef MPI
if (  sys.rank == 0 ){
  sys.vx=(double *)malloc(sys.natoms*sizeof(double));
  sys.vy=(double *)malloc(sys.natoms*sizeof(double));
  sys.vz=(double *)malloc(sys.natoms*sizeof(double));
  sys.fx=(double *)malloc(sys.natoms*sizeof(double));
  sys.fy=(double *)malloc(sys.natoms*sizeof(double));
  sys.fz=(double *)malloc(sys.natoms*sizeof(double));
}else{
  sys.fx=(double *)malloc(sys.my_range*sizeof(double));
  sys.fy=(double *)malloc(sys.my_range*sizeof(double));
  sys.fz=(double *)malloc(sys.my_range*sizeof(double));
}

if (  sys.rank == 0 ){

      sys.ranges = (int*)malloc(sys.size*sizeof(int));
      sys.disp = (int*)malloc(sys.size*sizeof(int));
      sys.ranges[0]=sys.my_range;
      sys.disp[0]=0;
      printf("ranges\tdisp\n%d\t%d\n", sys.ranges[0], sys.disp[0]);
      for (i=1; i<sys.size; i++){
        sys.ranges[ i ] = sys.natoms /  sys.size;
        if ( i < (sys.natoms % sys.size) )
          ++sys.ranges[ i ];
        sys.disp[i] = sys.disp[i-1]+sys.ranges[i];
        printf("%d\t%d\n", sys.ranges[i], sys.disp[i]);
      }
#else
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.natoms*sizeof(double));
    sys.fy=(double *)malloc(sys.natoms*sizeof(double));
    sys.fz=(double *)malloc(sys.natoms*sizeof(double));


#endif
    /* read restart */
    fp=fopen(restfile,"r");
    if(fp) {
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
        }
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
        }
        fclose(fp);
        azzero(sys.fx, sys.natoms);
        azzero(sys.fy, sys.natoms);
        azzero(sys.fz, sys.natoms);
    } else {
        perror("cannot read restart file");
        return 3;
    }

    t_start=cclock();
    /* initialize forces and energies.*/
    sys.nfi=0;
#ifdef MPI
}
#endif
    force(&sys);
#ifdef MPI
if (  sys.rank == 0 ){
#endif
    ekin(&sys);

    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);
#ifdef MPI
}
#endif
    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {
#ifdef MPI
if (  sys.rank == 0 ){
#endif
        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
        velverlet1(&sys);
#ifdef MPI
}
#endif
        force(&sys);
#ifdef MPI
if (  sys.rank == 0 ){
#endif
        velverlet2(&sys);

        ekin(&sys);
#ifdef MPI
}
#endif
    }
    /**************************************************/
#ifdef MPI
if (  sys.rank == 0 ){
#endif
    t_stop=cclock();
    /* clean up: close files, free memory */
    printf("Simulation Done. Elapsed time: %9.4f secs\n", t_stop - t_start );
    fclose(erg);
    fclose(traj);

    free( sys.disp );
    free( sys.ranges );
#ifdef MPI
}
#endif
    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);
#ifdef MPI
MPI_Finalize();
#endif
    return 0;
}
