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
    int nprint1, nprint2;
    char restfile1[BLEN], trajfile1[BLEN], ergfile1[BLEN], line[BLEN];
    char restfile2[BLEN], trajfile2[BLEN], ergfile2[BLEN];
    FILE *fp;
    mdsys_t sys1, sys2;
    fp = fopen("test.dat","w");

    if (fp) {
      sys1.natoms = 172;
      fprintf(fp, "%d\t\t# natoms\n", sys1.natoms);

      sys1.mass = 27.043;
      fprintf(fp, "%f\t# mass in AMU\n", sys1.mass);

      sys1.epsilon = 0.7823;
      fprintf(fp, "%f\t# epsilon in kcal/mol\n", sys1.epsilon);

      sys1.sigma = 7.251;
      fprintf(fp, "%f\t# sigma in angstrom\n", sys1.sigma);

      sys1.rcut = 6.3;
      fprintf(fp, "%f\t# rcut in angstrom\n", sys1.rcut);

      sys1.box = 16.582;
      fprintf(fp, "%f\t# box length (in angstrom)\n", sys1.box);

      strcpy(restfile1, "rest_test.rest");
      fprintf(fp, "%s\t# restart\n", restfile1);

      strcpy(trajfile1, "traj_test.xyz");
      fprintf(fp, "%s\t# trajectory\n", trajfile1);

      strcpy(ergfile1, "erg_test.dat");
      fprintf(fp, "%s\t# energies\n", ergfile1);

      sys1.nsteps = 59281;
      fprintf(fp, "%d\t\t# nr MD steps\n", sys1.nsteps);

      sys1.dt = 5.0;
      fprintf(fp, "%f\t# MD time step (in fs)\n", sys1.dt);

      nprint1 = 200;
      fprintf(fp, "%d\t\t# output print frequency\n", nprint1);

      fclose(fp);
    } else {
        perror("cannot read restart file");
        return 1;
    }

    fp = fopen("test.dat","r");
    /* read input file */
    if (fp){
      if(get_a_line(fp,line)) return 2;
      sys2.natoms=atoi(line);
      if(get_a_line(fp,line)) return 2;
      sys2.mass=atof(line);
      if(get_a_line(fp,line)) return 2;
      sys2.epsilon=atof(line);
      if(get_a_line(fp,line)) return 2;
      sys2.sigma=atof(line);
      if(get_a_line(fp,line)) return 2;
      sys2.rcut=atof(line);
      if(get_a_line(fp,line)) return 2;
      sys2.box=atof(line);
      if(get_a_line(fp,restfile2)) return 2;
      if(get_a_line(fp,trajfile2)) return 2;
      if(get_a_line(fp,ergfile2)) return 2;
      if(get_a_line(fp,line)) return 2;
      sys2.nsteps=atoi(line);
      if(get_a_line(fp,line)) return 2;
      sys2.dt=atof(line);
      if(get_a_line(fp,line)) return 2;
      nprint2=atoi(line);
    } else {
        printf("cannot read from file");
        return 1;
    }

    if ( sys1.natoms != sys2.natoms ) return 3;
    if ( sys1.mass != sys2.mass ) return 4;
    if ( sys1.epsilon != sys2.epsilon ) return 5;
    if ( sys1.sigma != sys2.sigma ) return 6;
    if ( sys1.rcut != sys2.rcut ) return 7;
    if ( sys1.box != sys2.box ) return 8;
    if ( strcmp(restfile1, restfile2) ) return 9;
    if ( strcmp(trajfile1, trajfile2) ) return 10;
    if ( strcmp(ergfile1, ergfile2) ) return 11;
    if ( sys1.nsteps != sys2.nsteps ) return 12;
    if ( sys1.dt != sys2.dt ) return 13;
    if ( nprint1 != nprint2 ) return 14;

    return 0;
}
