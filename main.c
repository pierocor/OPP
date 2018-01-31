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

    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.natoms*sizeof(double));
    sys.fy=(double *)malloc(sys.natoms*sizeof(double));
    sys.fz=(double *)malloc(sys.natoms*sizeof(double));

    /* define here the number of cells per dimension, using the integer division
    each cell must have side larger than RCUT*/
    sys.lc = sys.box/sys.rcut;  /* AAA */
    
    /*once i define the number of cells per dimension, i calculate the side of each cell with the double division*/
    sys.rc = sys.box/sys.lc;    /* AAA */
    
    /*total number of cell is to the power 3*/
    sys.ncells = sys.lc*sys.lc*sys.lc; /* AAA */
    printf("System divided in %d cells, of %f lenght.\n",sys.ncells, sys.rc);  /* AAA */
    
    /*each cell must be larger than RCUT, so if the total box L is not larger the 3*RCUT, or the total numbe rof cell is < 27, i don't use the cells
     i.e. the system is too small to benefit */
    sys.yescell = 1; /*integer flag to distinguish between the two cases*/  /* AAA */
    if(sys.ncells < 27) sys.yescell = 0; /* AAA */
    
    /*allocate the vector lscl which cointains the index of the particles,
     so that, the value stored in this array indicates the position in the same array of the next particle within the subcell */
    sys.lscl=(int *)malloc(sys.natoms*sizeof(int)); /* AAA */
    /*vector of dimension ncells which contains the index of HEAD particle in the array lscl */
    sys.head=(int *)malloc((sys.ncells)*sizeof(int)); /* AAA */

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
    /* refold all the positions in the box, i.e. coordinates must always be x,y,z in (0,L)*/
    if(sys.yescell) Putinthebox(&sys);  /*AAA*/
    force(&sys);
    ekin(&sys);

    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
        velverlet1(&sys);
       /* refold all the positions in the box, i.e. coordinates must always be x,y,z in (0,L)*/
       if(sys.yescell) Putinthebox(&sys); /* AAA */
        force(&sys);
        velverlet2(&sys);

        ekin(&sys);
    }
    /**************************************************/
    t_stop=cclock();
    /* clean up: close files, free memory */

    printf("Simulation Done. Elapsed time: %9.6f secs\n", t_stop - t_start );

    fclose(erg);
    fclose(traj);

    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);

    free(sys.lscl); /*AAA*/
    free(sys.head); /*AAA*/
    return 0;
}
