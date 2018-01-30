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
    int  i;
    FILE *fp;
    mdsys_t sys;

    /* pre-filled sys with selected values */
    sys.natoms = 3;
    sys.mass = 39.948;
    sys.epsilon = 0.2379;
    sys.sigma = 3.405;
    sys.rcut = 8.5;
    sys.box = 17.158;
    sys.nsteps = 1;
    sys.dt  = 5.0;


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

    /* read restart */
    fp=fopen("argon_time_step.rest","r");
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

    /* initialize forces and energies.*/
    sys.nfi=0;
    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("\tVx \t\tVy \t\tVz \n");
    for (i=0; i<sys.natoms; ++i) {
      printf("\t%f \t%f \t%f \n", sys.vx[i], sys.vy[i], sys.vz[i]);
      fprintf(fp, "\t%f \t%f \t%f \n", sys.vx[i], sys.vy[i], sys.vz[i]);
    }


    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

      /* propagate system and recompute energies */
      velverlet1(&sys);
      velverlet2(&sys);
      ekin(&sys);
    }
    /**************************************************/

    fp=fopen("single_time_step.dat","w");
    for (i=0; i<sys.natoms; ++i) {
      // printf("\t%f \t%f \t%f \n", sys.vx[i], sys.vy[i], sys.vz[i]);
      fprintf(fp, "\t%f \t%f \t%f \n", sys.vx[i], sys.vy[i], sys.vz[i]);
    }
    fclose(fp);

    /* clean up: close files, free memory */
    printf("Simulation Done.\n");

    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);

    return 0;
}
