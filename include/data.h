#ifndef DATA_H
#define DATA_H

/* generic file- or pathname buffer length */
#ifndef BLEN
#define BLEN 200
#endif

#ifndef EMPTY
#define EMPTY -1  /* AAA */
#endif


/* a few physical constants */
static const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
static const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;   
    int *lscl, *head;   /*AAA */
    int lc,ncells;   /*AAA */
    double rc;      /*AAA */
    int yescell;    /*AAA*/

};
typedef struct _mdsys mdsys_t;

#endif
