
#include <math.h>
#include <data.h>
#include <velocity_verlet.h>
// #include <../src/output.c>
// #include <../src/input.c>
// #include <../src/utilities.c>
// #include <../src/verlet1.c>
// #include <../src/verlet2.c>

/* compute forces */

__global__ void focres(double * d_rx,double *d_ry,double * d_rz, double * d_vx,double * d_vy,double * d_vz,double * d_fx,\
  double * d_fy,double * d_fz,double epot, double  epsilon,double sigma,int natoms, double box, double rcut){

    double ffac,r2,r6,r12,sigma2;
    double rx,ry,rz;
    // int i,j,n;

    sigma2 = sigma*sigma;

    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    int idy = threadIdx.y + blockIdx.y*blockDim.y;

    epot = 0.0;
     d_fx[idx] = 0.0;
     d_fy[idx] = 0.0;
     d_fz[idx] = 0.0;


        // if (i==j) continue;
        if (idx > idy){
            while ((d_rx[idx] - d_rx[idy]) >  0.5*box) rx -= box; while ((d_rx[idx] - d_rx[idy]) <  0.5*box) rx += box;
            while ((d_ry[idx] - d_ry[idy]) >  0.5*box) ry -= box; while ((d_ry[idx] - d_ry[idy]) <  0.5*box) ry += box;
            while ((d_rz[idx] - d_rz[idy]) >  0.5*box) rz -= box; while ((d_rz[idx] - d_rz[idy]) <  0.5*box) rz += box;

            r2 = rx*rx + ry*ry + rz*rz;

            if (r2 < rcut*rcut ) {
                r2 = (sigma2)/r2;
                r6 = r2*r2*r2;
                r12 = r6*r6;

                ffac = -4.0*epsilon*(-12.0*r12+6*r6)*r2/sigma2;
                epot += 4.0*epsilon*(r12-r6);

                d_fx[idx] += rx*ffac; d_fx[idy] -= rx*ffac;
                d_fy[idx] += ry*ffac; d_fy[idy] -= ry*ffac;
                d_fz[idx] += rz*ffac; d_fz[idy] -= rz*ffac;
            }
        }
    }
