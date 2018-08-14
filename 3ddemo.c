/* 3D simulation with dipole source at center of grid. */

#include "fdtd-alloc.h"
#include "fdtd-macro.h" 
#include "fdtd-proto.h"
#include "ezinc.h"
#include <stdio.h>
#include <complex.h>
void runsim(double k_x, double k_y, double k_z, double t_seed, double i)
{
  Grid *g;
  
  ALLOC_1D(g, 1, Grid); // allocate memory for grid structure
  MaxTime = t_seed;
  double Kx = k_x;
  double Ky = k_y;
  double Kz = k_z;
  double delta1, delta2;
  printf("Kx is %f,  Ky is %f, Kz is %f \n", Kx, Ky, Kz);
  int j = (int) i;
  sensorInit(g, Kx, Ky, Kz, j);
  gridInit(g, Kx, Ky, Kz);        // initialize 3D grid
  //abcInit(g);         // initialize ABC
  ezIncInit(g); 
  //snapshot3dInit(g);  // initialize snapshots
  /* do time stepping */
   
  for (Time = 0; Time < MaxTime; Time++) {
    updateH(g);       // update magnetic fields 
    ezInc(g, Time, 5, 13,  1, Kx, Ky, Kz);//*cexp(I*(Kx*16 + Ky*16)*(1e-8));
    //   ezInc(g, Time, 36, 11,  1, Kx, Ky, Kz);
    //ezInc(g, Time, 12, 30,  1, Kx, Ky, Kz);

    
    //*cexp(I*(Kx*66 + Ky*103)*(1e-8));
    updateE(g);       // update electric fields 
    // abc(g);           // apply ABC
    // snapshot3d(g);    // take a snapshot (if appropriate)
    Transmission(g, j, Time);
  } // end of time-stepping


}
