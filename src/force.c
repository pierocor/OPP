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
void force(mdsys_t *sys){
    double ffac,r2,r6,r12,sigma2;
    double rx,ry,rz;

    int i,j;
#ifdef CL
    int a,lcyz,lcxyz,mc[3],c,mc1[3],c1;
    double rshift[3];
#endif
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
#else
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
#endif
#ifdef CL

    if(sys->yescell) {

    lcyz = sys->lc*sys->lc;
    lcxyz = sys->lc*lcyz;

    /* Reset the  head */
    for (c=0; c<lcxyz; c++) sys->head[c] = EMPTY;

    /* Scan atoms to construct head and lscl  */
    for (i=0; i<(sys->natoms); i++) {
        mc[0] = sys->rx[i]/sys->rc;
        mc[1] = sys->ry[i]/sys->rc;
        mc[2] = sys->rz[i]/sys->rc;
        /*here is essential that the coordinates rx,ry,rz are in the box (0,L)
         otherwise i can not assing them to the correct cell, using mc[] */

        //if(sys->rx[i]/sys->rc > sys->lc) printf("X %f %f %d \n", sys->rx[i], sys->rx[i]/sys->rc , mc[0] );
        //if(sys->ry[i]/sys->rc > sys->lc) printf("Y %f %f %d \n", sys->ry[i], sys->ry[i]/sys->rc , mc[1] );
        //if(sys->rz[i]/sys->rc > sys->lc) printf("Z %f %f %d \n", sys->rz[i], sys->rz[i]/sys->rc , mc[2] );

        /* Translate the vector cell index, mc, to a scalar cell index AAAA*/
        c = mc[0]*lcyz + mc[1]*sys->lc + mc[2];

        /* Link to the previous occupant (or EMPTY if you're the 1st) */
        sys->lscl[i] = sys->head[c];

        /* The last one goes to the header  */
        sys->head[c] = i;
    } /* Endfor atom i  */

    /* Scan inner cells */
    for (mc[0]=0; mc[0]<sys->lc; (mc[0])++)
    for (mc[1]=0; mc[1]<sys->lc; (mc[1])++)
    for (mc[2]=0; mc[2]<sys->lc; (mc[2])++) {

                /* Calculate a scalar cell index */
                c = mc[0]*lcyz+mc[1]*sys->lc+mc[2];
                /* Skip this cell if empty */
                if (sys->head[c] == EMPTY) continue;

                /* Scan the neighbor cells (including itself) of cell c */
                for (mc1[0]=mc[0]-1; mc1[0]<=mc[0]+1; (mc1[0])++)
                for (mc1[1]=mc[1]-1; mc1[1]<=mc[1]+1; (mc1[1])++)
                for (mc1[2]=mc[2]-1; mc1[2]<=mc[2]+1; (mc1[2])++) {

                           /*if the nearest neighbour cell index is outside the box, i remember this and use
                            shifted cordinates (by -L or + L) later where I evaluate particle distances */
                           for (a=0; a<3; a++) {
                               if (mc1[a] < 0)
                                rshift[a] = -sys->box;
                              else if (mc1[a]>=sys->lc)
                                rshift[a] = sys->box;
                              else
                                rshift[a] = 0.0;
                            }
                            /* Calculate the scalar cell index of the neighbor cell */
                            /* by using the modulus(%) the cell is correcly identified with PBC,
                             but the position of the particles inside need to be trnaslated usign rshift[] */
                            c1 = ((mc1[0]+sys->lc)%sys->lc)*lcyz+((mc1[1]+sys->lc)%sys->lc)*sys->lc+((mc1[2]+sys->lc)%sys->lc);

                            /* Skip this neighbor cell if empty */
                            if (sys->head[c1] == EMPTY) continue;

                            /* Scan atom i in cell c */
                            i = sys->head[c];
                            while (i != EMPTY) {

                                /* Scan atom j in cell c1 */
                                j = sys->head[c1];
                                while (j != EMPTY) {

                                    /* Avoid double counting of pairs */
                                    if (i < j) {
                                        /* here i dont use the periodic image because the subcell take care of this using the rshift[]
                                         defined above*/
                                        rx=sys->rx[i] - ( sys->rx[j]+rshift[0]);
                                        ry=sys->ry[i] - ( sys->ry[j]+rshift[1]);
                                        rz=sys->rz[i] - ( sys->rz[j]+rshift[2]);
                                        r2 = rx*rx+ry*ry+rz*rz;
                                        if (r2 < sys->rcut*sys->rcut) {

                                            r2= (sigma2)/r2;
                                            r6= r2*r2*r2;
                                            r12= r6*r6;


                                            ffac = -4.0*sys->epsilon*(-12.0*r12+6*r6)*r2 /sigma2;

                                            sys->epot += 4.0*sys->epsilon*(r12-r6);


                                            sys->fx[i] += rx*ffac; sys->fx[j] -= rx*ffac;
                                            sys->fy[i] += ry*ffac; sys->fy[j] -= ry*ffac;
                                            sys->fz[i] += rz*ffac; sys->fz[j] -= rz*ffac;


                                            }
                                    } /* Endif i<j */

                                    j = sys->lscl[j];
                                } /* Endwhile j not empty */

                                i = sys->lscl[i];
                            } /* Endwhile i not empty */

                        } /* Endfor neighbor cells, c1 */
            } /* Endfor central cell, c */


    }
    else{  /*if the system is so small i use the old code without cells*/

        /*old code */
        for(i=0; i < (sys->natoms); ++i) {
            for(j=i+1; j < (sys->natoms); ++j) {

                if (i==j) continue;

                /* get distance between particle i and j */
                rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
                ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
                rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
                r2 = rx*rx+ry*ry+rz*rz;

                if (r2 < sys->rcut*sys->rcut) {


                    r2= (sigma2)/r2;
                    r6= r2*r2*r2;
                    r12= r6*r6;

                    ffac = -4.0*sys->epsilon*(-12.0*r12+6*r6)*r2/sigma2;

                    sys->epot += 4.0*sys->epsilon*(r12-r6);

                    sys->fx[i] += rx*ffac; sys->fx[j] -= rx*ffac;
                    sys->fy[i] += ry*ffac; sys->fy[j] -= ry*ffac;
                    sys->fz[i] += rz*ffac; sys->fz[j] -= rz*ffac;
                }
            }
        }



    } //yescell
}
#else  // IFNDEF CL
#ifdef MPI
    #ifdef OPENMP
    double epot = 0.0;
    #pragma omp parallel for default(shared) schedule(guided) private(ii, i, ffac, r2, r6, r12, rx, ry, rz, j) reduction(+:epot)
    #endif
    for(ii = 0; ii < (sys->natoms); ii += sys->size) {
      i = ii + sys->rank;
      if ( i  >= (sys->natoms ) )
        continue;
#else
    #ifdef OPENMP
    double epot = 0.0;
    #pragma omp parallel for default(shared) schedule(guided) private(i, ffac, r2, r6, r12, rx, ry, rz, j) reduction(+:epot)
    #endif
    for(i=0; i < (sys->natoms); ++i) {
#endif
      for(j=0; j < (sys->natoms); ++j) {

        if ( i == j ) continue;
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
          epot += 0.5 * 4.0*sys->epsilon*(r12-r6); /*REDEFINED epot, NO MORE pow */
          #else
          sys->epot += 0.5 * 4.0 *sys->epsilon*(r12-r6); /*REDEFINED epot, NO MORE pow */
          #endif
#ifdef MPI
          sys->cx[i] += rx*ffac;
          sys->cy[i] += ry*ffac;
          sys->cz[i] += rz*ffac;
#else
          sys->fx[i] += rx*ffac;
          sys->fy[i] += ry*ffac;
          sys->fz[i] += rz*ffac;
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
#endif // ifdef CL
