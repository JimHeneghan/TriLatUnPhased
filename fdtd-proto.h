#ifndef _FDTD_PROTO_H
#define _FDTD_PROTO_H
#include <complex.h>
#include "fdtd-grid.h"

/* Function prototypes */
void abcInit(Grid *g);
void abc(Grid *g);

void gridInit(Grid *g, double Kx, double Ky, double Kz);

void snapshot3dInit(Grid *g);
void snapshot3d(Grid *g);

void updateE(Grid *g);
void updateH(Grid *g);


void sensorInit(Grid *g, double Kx, double Ky, double Kz, int j);
void Transmission(Grid *g, int j,  double time);
#endif
