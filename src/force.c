
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
    double r,ffac;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
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
      for(j=0; j < (sys->natoms); ++j) {

        /* particles have no interactions with themselves */
        if (i==j) continue;

        /* get distance between particle i and j */
        rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
        ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
        rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
        r = sqrt(rx*rx + ry*ry + rz*rz);

        /* compute force and energy if within cutoff */
        if (r < sys->rcut) {
          ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                   +6*pow(sys->sigma/r,6.0)/r);

          sys->epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                         -pow(sys->sigma/r,6.0));
#ifdef MPI
          sys->cx[i] += rx/r*ffac;
          sys->cy[i] += ry/r*ffac;
          sys->cz[i] += rz/r*ffac;
#else
          sys->fx[i] += rx/r*ffac;
          sys->fy[i] += ry/r*ffac;
          sys->fz[i] += rz/r*ffac;
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
