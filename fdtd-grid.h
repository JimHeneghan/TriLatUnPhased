#ifndef _FDTD_GRID_H
#define _FDTD_GRID_H
#include <complex.h>
enum GRIDTYPE {oneDGrid, teZGrid, tmZGrid, threeDGrid};

struct Grid {
  complex double *hx, *hy, *hz;
  double *chxh, *chxe;
  double *chyh, *chye;
  double *chzh, *chze;
  complex double *ex, *ey, *ez;
  double *cexe, *cexh;
  double *ceye, *ceyh;
  double *ceze, *cezh;
  double *epsR;
  int sizeX, sizeY, sizeZ;
  int time, maxTime;
  int type;
  double cdtds;
  complex double phix, phiy, phiz;
};

typedef struct Grid Grid;

#endif
