
#include <math.h>
#include <data.h>
#include <velocity_verlet.h>
#ifdef MPI
#include <mpi.h>
#endif

#ifdef OPENMP
#include <omp.h>
#endif

/* helper function: apply minimum image convention */
static double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

/* compute forces */

void force(mdsys_t *sys)
{
    double ffac,r2,r6,r12,sigma2;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    sigma2=sys->sigma*sys->sigma;  /* <-DEFINED NEW HERE OUTSIDE LOOP */
#ifdef MPI
    int ii;

    /* BROADCAST OF POSITION AND VELOCITY VECTORS*/
    MPI_Bcast( sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast( sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast( sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    azzero(sys->cx,sys->natoms);
    azzero(sys->cy,sys->natoms);
    azzero(sys->cz,sys->natoms);
    #ifdef OPENMP
    double epot = 0.0;
    #pragma omp parallel for default(shared) private(ii, i, ffac, r2, r6, r12, rx, ry, rz, j) reduction(+:epot)
    #endif
    for(ii = 0; ii < (sys->natoms - 1); ii += sys->size) {
      i = ii + sys->rank;
      if ( i  >= (sys->natoms - 1) )
        continue;
#else
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
    #ifdef OPENMP
    double epot = 0.0;
    #pragma omp parallel for default(shared) private(i, ffac, r2, r6, r12, rx, ry, rz, j) reduction(+:epot)
    #endif
    for(i=0; i < (sys->natoms); ++i) {
#endif
      for(j=i+1; j < (sys->natoms); ++j) {
        /* get distance between particle i and j */
        rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
        ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
        rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
        r2 = rx*rx + ry*ry + rz*rz;    /* <-DEFINED R2 HERE TO DELETE SQRT instead of r = sqrt(rx*rx + ry*ry + rz*rz)*/

        if (r2 < sys->rcut*sys->rcut ) {   /* changed if (r < sys->rcut) */

          r2 = (sigma2)/r2;  /* <-REDEFINED R2 HERE */
          r6 = r2*r2*r2; /* <-DEFINED R6 HERE */
          r12 = r6*r6; /* <-DEFINED R12 HERE */

          ffac = -4.0*sys->epsilon*(-12.0*r12+6*r6)*r2/sigma2; /*<-REDEFINED ffac,NO MORE pow*/
          #ifdef OPENMP
          epot += 4.0*sys->epsilon*(r12-r6); /*REDEFINED epot, NO MORE pow */
          #else
          sys->epot += 4.0*sys->epsilon*(r12-r6); /*REDEFINED epot, NO MORE pow */
          #endif
#ifdef MPI
          sys->cx[i] += rx*ffac; sys->cx[j] -= rx*ffac;
          sys->cy[i] += ry*ffac; sys->cy[j] -= ry*ffac;
          sys->cz[i] += rz*ffac; sys->cz[j] -= rz*ffac;
#else
          #ifdef OPENMP
          #pragma omp atomic update
          #endif
          sys->fx[i] += rx*ffac;
          #ifdef OPENMP
          #pragma omp atomic update
          #endif
          sys->fx[j] -= rx*ffac;
          #ifdef OPENMP
          #pragma omp atomic update
          #endif
          sys->fy[i] += ry*ffac;
          #ifdef OPENMP
          #pragma omp atomic update
          #endif
          sys->fy[j] -= ry*ffac;
          #ifdef OPENMP
          #pragma omp atomic update
          #endif
          sys->fz[i] += rz*ffac;
          #ifdef OPENMP
          #pragma omp atomic update
          #endif
          sys->fz[j] -= rz*ffac;
#endif
        }
      }
    }
    #ifdef OPENMP
    sys->epot = epot;
    #endif
#ifdef MPI

    MPI_Reduce( sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    double tmp = sys->epot;
    MPI_Reduce ( &tmp, &(sys->epot), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
}
