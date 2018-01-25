/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <data.h>
#include <velocity_verlet.h>



/* a few physical constants */
extern const double kboltz;
extern const double mvsq2e;

/* main */
int main(int argc, char **argv)
{
    int i;
    // char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp;
    mdsys_t sys;

    /* pre-filled sys with selected values */
    sys.natoms = 2916;
    sys.mass = 39.948;
    sys.epsilon = .2379;
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


    /* read restart */
    fp=fopen("argon_2916.rest","r");
    if(fp) {
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
        }
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
        }
        fclose(fp);
    } else {
        perror("cannot read restart file");
        return 3;
    }

    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);

    ekin(&sys);

    printf("Kinetic Energy of the system :\n");
    fp=fopen("kinetic_test.dat","w");
    printf("%f\n", sys.ekin);
    fprintf(fp, "%f\n", sys.ekin);
    fclose(fp);

    /* clean up: close files, free memory */
    printf("Simulation Done.\n");
    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);


    return 0;
}
