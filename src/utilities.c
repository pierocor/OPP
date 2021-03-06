
#include <stdio.h>
#include <stdlib.h>
#include <data.h>
#include <velocity_verlet.h>
#include <sys/time.h>
/* helper function: zero out an array */
void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

/* helper function: apply minimum image convention
double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
<<<<<<< HEAD
}*/

/* compute kinetic energy */
void ekin(mdsys_t *sys)
{
    int i;

    sys->ekin=0.0;
    for (i=0; i<sys->natoms; ++i) {
        sys->ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}

double cclock()
  /* Returns elepsed seconds past from the last call to timer rest */
{

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}


double SignR(double v,double x) {if (x > 0) return v; else return -v;} /*user defined sign fucntion */ /*AAA*/

/*fucntion to put particles in the box corectly
 if the one position is x > L, the x = X-L*/ /*AAA*/
void Putinthebox(mdsys_t *sys) {         /*AAA*/
    int i;                                /*AAA*/
    for (i=0; i<(sys->natoms); i++) {  /*AAA*/
    sys->rx[i] = sys->rx[i] - SignR(0.5*sys->box,sys->rx[i]) - SignR(0.5*sys->box,sys->rx[i]-sys->box);  /*AAA*/
    sys->ry[i] = sys->ry[i] - SignR(0.5*sys->box,sys->ry[i]) - SignR(0.5*sys->box,sys->ry[i]-sys->box);   /*AAA*/
    sys->rz[i] = sys->rz[i] - SignR(0.5*sys->box,sys->rz[i]) - SignR(0.5*sys->box,sys->rz[i]-sys->box);    /*AAA*/
    }
}      /*AAA*/
