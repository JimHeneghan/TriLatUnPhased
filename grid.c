#include "fdtd-macro.h"
#include "fdtd-alloc.h"
#include <complex.h>
#include <math.h>

void gridInit(Grid *g, double Kx, double Ky, double Kz) {
  double imp0 = 377.0, dx; 
  int mm, nn, pp;
  double XCenter1, XCenter2, YCenter1, YCenter2, r2, XLocC, YLocC, rad, XLoc1, XLoc2, YLoc1, YLoc2, dist1, dist2;
  double epsr;
  Type = threeDGrid;   /*@ \label{grid3dhomoA} @*/
  SizeX = 87; // size of domain
  SizeY = 50;
  SizeZ = 3;
  Cdtds = 1.0/2; // Courant number /*@ \label{grid3dhomoB} @*/
  dx = 0.01;
  Phix = cexp(I*Kx*SizeX*dx);
  Phiy = cexp(I*Ky*SizeY*dx);
  Phiz = cexp(I*Kz*SizeY*dx);

  /* memory allocation */
  ALLOC_3D(g->hx,   SizeX, SizeY, SizeZ, complex double); /*@ \label{grid3dhomoC} @*/
  ALLOC_3D(g->chxh, SizeX, SizeY, SizeZ, double);
  ALLOC_3D(g->chxe, SizeX, SizeY, SizeZ, double);
  ALLOC_3D(g->hy,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->chyh, SizeX, SizeY, SizeZ, double);
  ALLOC_3D(g->chye, SizeX, SizeY, SizeZ, double);
  ALLOC_3D(g->hz,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->chzh, SizeX, SizeY, SizeZ, double);
  ALLOC_3D(g->chze, SizeX, SizeY, SizeZ, double);

  ALLOC_3D(g->ex,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->cexe, SizeX, SizeY, SizeZ, double);
  ALLOC_3D(g->cexh, SizeX, SizeY, SizeZ, double);
  ALLOC_3D(g->ey,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->ceye, SizeX, SizeY, SizeZ, double);
  ALLOC_3D(g->ceyh, SizeX, SizeY, SizeZ, double);
  ALLOC_3D(g->ez,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->ceze, SizeX, SizeY, SizeZ, double);
  ALLOC_3D(g->cezh, SizeX, SizeY, SizeZ, double); /*@ \label{grid3dhomoD} @*/

  ALLOC_3D(g->epsR, SizeX, SizeY, SizeZ,  double);
  
  rad = 12; 
  XCenter1 = 22;
  YCenter1 = 12;
  XCenter2 = 65;
  YCenter2 = 37;
  pp = 1;
  r2 = rad*rad;
  /* for (mm = 0; mm < SizeX; mm++){ */
  /* 	XLoc1 = mm - XCenter1; */
  /*       XLoc2 = mm - XCenter2; */
  /* 	//       XLocC = mm - XCenter; */
  /*     for (nn = 0; nn < SizeY; nn++){  */
  /*       YLoc1 = nn - YCenter1; */
  /*       YLoc2 = nn - YCenter2; */
  /* 	//  YLocC = nn - YCenter; */
  /* 	for (pp = 0; pp < SizeZ; pp++){ */
  /* 	  if(XLoc1*XLoc1 + YLoc1*YLoc1 < r2 ) */
  /* 	    EpsR(mm, nn, pp) = 14; */
  /* 	  else if(XLoc2*XLoc2 + YLoc2*YLoc2 < r2) */
  /* 	    EpsR(mm, nn, pp) = 14; */
  /* 	  else */
  /* 	    EpsR(mm, nn, pp) = 1.0; */
  /* 	} */

  /*     } */
  /* } */
  epsr = 14;

  /*Zeroing fields*/

    for (mm = 0; mm < SizeX; mm++)
      for (nn = 0; nn < SizeY; nn++) 
	for (pp = 0; pp < SizeZ; pp++) {
	  Ex(mm, nn,pp) = 0;
	  Ey(mm, nn,pp) = 0;
	  Ez(mm, nn,pp) = 0;

	  Hx(mm, nn,pp) = 0;
	  Hy(mm, nn,pp) = 0;
	  Hz(mm, nn,pp) = 0;
	}
  /* set electric-field update coefficients */
  for (mm = 0; mm < SizeX; mm++) /*@ \label{grid3dhomoE} @*/
    for (nn = 0; nn < SizeY; nn++) 
      for (pp = 0; pp < SizeZ; pp++) {
	Cexe(mm, nn, pp) = 1.0;
	Cexh(mm, nn, pp) = Cdtds * imp0;
      }

  for (mm = 0; mm < SizeX; mm++)
    for (nn = 0; nn < SizeY; nn++) 
      for (pp = 0; pp < SizeZ; pp++) {
	Ceye(mm, nn, pp) = 1.0;
	Ceyh(mm, nn, pp) = Cdtds * imp0;
      }

  for (mm = 0; mm < SizeX; mm++)
    for (nn = 0; nn < SizeY; nn++) 
      for (pp = 0; pp < SizeZ; pp++) {
	Ceze(mm, nn, pp) = 1.0;
	Cezh(mm, nn, pp) = Cdtds * imp0;
      }

  /* set cylinders*/

   for (mm = 0; mm < SizeX; mm++) /*@ \label{grid3dhomoE} @*/
    for (nn = 0; nn < SizeY; nn++) 
      for (pp = 0; pp < SizeZ; pp++) {

	XLoc1 = XCenter1 - mm; 
        XLoc2 = XCenter2 - mm;
	YLoc1 = YCenter1 - nn;
	YLoc2 = YCenter2 - nn;
	dist1 = (int)(sqrt(XLoc1*XLoc1 + YLoc1*YLoc1));
	dist2 = (int)(sqrt(XLoc2*XLoc2 + YLoc2*YLoc2));
	if((dist1 < rad)||(dist2 < rad)){
	  Cexe(mm, nn, pp) = 1.0;
	  Cexh(mm, nn, pp) = Cdtds * imp0/epsr;
	}
      }

  for (mm = 0; mm < SizeX; mm++)
    for (nn = 0; nn < SizeY; nn++) 
      for (pp = 0; pp < SizeZ; pp++) {
	XLoc1 = XCenter1 - mm; 
        XLoc2 = XCenter2 - mm;
	YLoc1 = YCenter1 - nn;
	YLoc2 = YCenter2 - nn;
	dist1 = (int)(sqrt(XLoc1*XLoc1 + YLoc1*YLoc1));
	dist2 = (int)(sqrt(XLoc2*XLoc2 + YLoc2*YLoc2));
	if((dist1 < rad)||(dist2 < rad)){
	  Ceye(mm, nn, pp) = 1.0;
	  Ceyh(mm, nn, pp) = Cdtds * imp0/epsr;
	}
      }

  for (mm = 0; mm < SizeX; mm++)
    for (nn = 0; nn < SizeY; nn++) 
      for (pp = 0; pp < SizeZ; pp++) {
	XLoc1 = XCenter1 - mm; 
        XLoc2 = XCenter2 - mm;
	YLoc1 = YCenter1 - nn;
	YLoc2 = YCenter2 - nn;
	dist1 = (int)(sqrt(XLoc1*XLoc1 + YLoc1*YLoc1));
	dist2 = (int)(sqrt(XLoc2*XLoc2 + YLoc2*YLoc2));
	if((dist1 < rad)||(dist2 < rad)){
	  Ceze(mm, nn, pp) = 1.0;
	  Cezh(mm, nn, pp) = Cdtds * imp0/epsr;
	}
      }
  
  /* set magnetic-field update coefficients */
  for (mm = 0; mm < SizeX; mm++)
    for (nn = 0; nn < SizeY; nn++) 
      for (pp = 0; pp < SizeZ; pp++) {
	Chxh(mm, nn, pp) = 1.0;
	Chxe(mm, nn, pp) = Cdtds / imp0;
      }

  for (mm = 0; mm < SizeX; mm++)
    for (nn = 0; nn < SizeY; nn++) 
      for (pp = 0; pp < SizeZ; pp++) {
	Chyh(mm, nn, pp) = 1.0;
	Chye(mm, nn, pp) = Cdtds / imp0;
      }

  for (mm = 0; mm < SizeX; mm++)
    for (nn = 0; nn < SizeY; nn++) 
      for (pp = 0; pp < SizeZ; pp++) {
	Chzh(mm, nn, pp) = 1.0;
	Chze(mm, nn, pp) = Cdtds / imp0;
      } /*@ \label{grid3dhomoF} @*/

    float temp;
  char filename[100];
  FILE *out;
  sprintf(filename,"Media/Square.dat");
  out = fopen(filename,"w");
  pp = 1;
  for (nn=SizeY-1; nn>=0; nn--)
    for (mm=0; mm<SizeX; mm++) {
      temp = (float)Cezh(mm,nn, pp); // store data as a float
      fprintf(out, "%f \n", temp);
    }

  fclose(out);

  

  return;
}  /* end gridInit() */
