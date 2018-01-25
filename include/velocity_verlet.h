#ifndef VELOCITY_VERLET_H
#define VELOCITY_VERLET_H

/*the list of functions.*/
static void velverlet1(mdsys_t *sys)
static void velverlet2(mdsys_t *sys)
static void force(mdsys_t *sys)
static int get_a_line(FILE *fp, char *buf)
static void output(mdsys_t *sys, FILE *erg, FILE *traj)
static void azzero(double *d, const int n)
static double pbc(double x, const double boxby2)
static void ekin(mdsys_t *sys)

#endif
