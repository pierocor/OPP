
#include <math.h>
#include <data.h>
#include <velocity_verlet.h>



/* compute forces */

void force(mdsys_t *sys) 
{
    double ffac,r2,r6,r12,sigma2;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
    
    sigma2=sys->sigma*sys->sigma;  /* <-DEFINED NEW HERE OUTSIDE LOOP */

    for(i=0; i < (sys->natoms); ++i) {
        for(j=i+1; j < (sys->natoms); ++j) {

            if (i==j) continue;
            
            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            r2 = rx*rx + ry*ry + rz*rz;    /* <-DEFINED R2 HERE TO DELETE SQRT instead of r = sqrt(rx*rx + ry*ry + rz*rz)*/
            
            if (r2 < sys->rcut*sys->rcut ) {   /* CHANGED if (r < sys->rcut) */

            
                r2= (sigma2)/r2;  /* <-REDEFINED R2 HERE */
                r6= r2*r2*r2; /* <-DEFINED R6 HERE */
                r12= r6*r6; /* <-DEFINED R12 HERE */
            
            
                ffac = -4.0*sys->epsilon*(-12.0*r12+6*r6)*r2/sigma2; /* <-REDEFINED ffac, NO MORE pow*/
                
                sys->epot += 4.0*sys->epsilon*(r12-r6); /*REDEFINED epot, NO MORE pow */


                sys->fx[i] += rx*ffac; sys->fx[j] -= rx*ffac;
                sys->fy[i] += ry*ffac; sys->fy[j] -= ry*ffac;
                sys->fz[i] += rz*ffac; sys->fz[j] -= rz*ffac;
            }
        }
    }
}
