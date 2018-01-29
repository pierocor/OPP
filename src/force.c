
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

    if ( sys->rank == 0 ){
      azzero(sys->fx, sys->natoms);
      azzero(sys->fy, sys->natoms);
      azzero(sys->fz, sys->natoms);
    }else{
      azzero(sys->fx, sys->my_range);
      azzero(sys->fy, sys->my_range);
      azzero(sys->fz, sys->my_range);
    }
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
          sys->fx[i - sys->my_start] += rx/r*ffac;
          sys->fy[i - sys->my_start] += ry/r*ffac;
          sys->fz[i - sys->my_start] += rz/r*ffac;
#else
          sys->fx[i] += rx/r*ffac;
          sys->fy[i] += ry/r*ffac;
          sys->fz[i] += rz/r*ffac;
#endif
        }
      }
    }
#ifdef MPI
    // if( sys->rank != 0 ){
    //   MPI_Send(sys->fx, sys->my_range, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD);
    //   MPI_Send(sys->fy, sys->my_range, MPI_DOUBLE, 0, TAG+1, MPI_COMM_WORLD);
    //   MPI_Send(sys->fz, sys->my_range, MPI_DOUBLE, 0, TAG+2, MPI_COMM_WORLD);
    // }else{
    //   int section_size = sys->my_range;
    //   int position = 0;
    //   for ( int sender = 1; sender <  sys->size; sender++ ){
    //     position += section_size;
    //     if ( sender == sys->natoms %  sys->size )
    //       section_size--;
    //     MPI_Recv( &(sys->fx[ position ]), section_size, MPI_DOUBLE, sender, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     MPI_Recv( &(sys->fy[ position ]), section_size, MPI_DOUBLE, sender, TAG+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     MPI_Recv( &(sys->fz[ position ]), section_size, MPI_DOUBLE, sender, TAG+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //   }
    // }
    ////////////////////////
    if( sys->rank != 0 ){
      MPI_Gatherv( sys->fx, sys->my_range, MPI_DOUBLE, sys->fx, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv( sys->fy, sys->my_range, MPI_DOUBLE, sys->fy, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv( sys->fz, sys->my_range, MPI_DOUBLE, sys->fz, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }else{
      MPI_Gatherv( sys->fx, sys->my_range, MPI_DOUBLE, sys->fx, sys->ranges, sys->disp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv( sys->fy, sys->my_range, MPI_DOUBLE, sys->fy, sys->ranges, sys->disp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv( sys->fz, sys->my_range, MPI_DOUBLE, sys->fz, sys->ranges, sys->disp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    double tmp = (sys->epot);
    MPI_Reduce ( &tmp, &(sys->epot), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
}
