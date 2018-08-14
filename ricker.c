#include "fdtd-alloc.h"
#include "fdtd-macro.h" 
#include "fdtd-proto.h"
#include "ezinc.h"
#include <stdio.h>
#include <complex.h>
#include <math.h>
static double cdtds, ppw = 0;

/* initialize source-function variables */
void ezIncInit(Grid *g){

  ppw = 30;
  cdtds = Cdtds;
  return;
}

/* calculate source function at given time and location */
void ezInc(Grid *g, int time, double LocX, double LocY, double LocZ, double Kx, double Ky, double Kz) {
  double arg, dx, je, omega, dt, delta, deltat, arg2;
  double delx, dely, delz;
  double LocX2, LocY2, LocZ2;
  delx = 43;
  dely = 25.0;
  delz = 0.0;
  double dTime = (double) time; 
  LocX2 = LocX + delx;
  LocY2 = LocY + dely;
  LocZ2 = LocZ + delz;
  dx = 0.01;
  je = 50.0;
  dt = dx/(3e8*2.0);
  omega  = 2.0*3e8*M_PI/(je*dx);

   
  arg = (dTime*dt - dt*1000.0)/50.0/dt;
  arg2 = arg * arg;
  delta =  (delx*Kx + dely*Ky + delz*Kz)*dx;

  Hz((int)LocX, (int)LocY, (int)LocZ) = Hz((int)LocX, (int)LocY, (int)LocZ) - exp(-0.5*arg2)*sin(omega*(dTime*dt - 1000*dt));
  Hz((int)LocX2, (int)LocY2, (int)LocZ2) = Hz((int)LocX2, (int)LocY2, (int)LocZ2) - exp(-0.5*arg2)*sin(omega*(dTime*dt - 1000*dt))*cexp(I*delta);

  
}
