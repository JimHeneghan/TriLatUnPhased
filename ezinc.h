#ifndef _EZINC_H
#define _EZINC_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "fdtd-macro.h"
#include <complex.h>

void ezIncInit(Grid *g); 
void ezInc(Grid *g, int time, double LocX, double LocY, double LocZ, double Kx, double Ky, double Kz);


#endif
