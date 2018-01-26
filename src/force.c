
#include <math.h>
#include <data.h>
#include <velocity_verlet.h>



/* compute forces */
void force(mdsys_t *sys) 
{
    double r,ffac;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    for(i=0; i < (sys->natoms); ++i) {
        for(j=i+1; j < (sys->natoms); ++j) {   /* <-CHANGED HERE j=1+1 instead of j=0 ...in this way I do not have to make a complete loop on i and then a complete loop on j counting twice the i <--> j interaction. Now I count it ones.*/ 

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
                
                sys->epot += 4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                               -pow(sys->sigma/r,6.0));  /* <-CHANGED HERE . I multiply by 2 the contribution of the potential energy epot because now i count each couple ones ..removed 0.5*4.0*sys... */

                sys->fx[i] += rx/r*ffac; sys->fx[j] -= rx/r*ffac;  /* <-CHANGED HERE added f[j] */
                sys->fy[i] += ry/r*ffac; sys->fy[j] -= ry/r*ffac;  /* <-CHANGED HERE added f[j] */
                sys->fz[i] += rz/r*ffac; sys->fz[j] -= rz/r*ffac;  /* <-CHANGED HERE added f[j] */


            }
        }
    }
}

