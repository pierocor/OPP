
#include <math.h>
#include <data.h>
#include <velocity_verlet.h>
#ifdef MPI
#include <mpi.h>
#define TAG 50
#endif


/* compute forces */

void force(mdsys_t *sys)
{
    double ffac,r2,r6,r12,sigma2;
    double rx,ry,rz;
    // int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    sigma2=sys->sigma*sys->sigma;  /* <-DEFINED NEW HERE OUTSIDE LOOP */
#ifdef MPI
    /* BROADCAST OF POSITION AND VELOCITY VECTORS*/
    MPI_Bcast( sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast( sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast( sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    azzero(sys->cx,sys->natoms);
    azzero(sys->cy,sys->natoms);
    azzero(sys->cz,sys->natoms);
    for(i = sys->my_start; i < sys->my_end; ++i ) {
#else
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
    for(i=0; i < (sys->natoms); ++i) {
#endif
      for(int j=i+1; j < (sys->natoms); ++j) {

        if (i==j) continue;

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

          sys->epot += 4.0*sys->epsilon*(r12-r6); /*REDEFINED epot, NO MORE pow */
#ifdef MPI
          sys->cx[i] += rx*ffac; sys->cx[j] -= rx*ffac;
          sys->cy[i] += ry*ffac; sys->cy[j] -= ry*ffac;
          sys->cz[i] += rz*ffac; sys->cz[j] -= rz*ffac;
#else
          sys->fx[i] += rx*ffac; sys->fx[j] -= rx*ffac;
          sys->fy[i] += ry*ffac; sys->fy[j] -= ry*ffac;
          sys->fz[i] += rz*ffac; sys->fz[j] -= rz*ffac;
#endif
        }
      }
    }
#ifdef MPI

    MPI_Reduce( sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    double tmp = sys->epot;
    MPI_Reduce ( &tmp, &(sys->epot), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
}
