#include <stdio.h>
#include <stdlib.h>
#include "fdtd-macro.h"
#include <complex.h>
#include <math.h>
static int temporalStride = -2, frameX = 0, frameY = 0, startTime;
static char basename[80];

void snapshot3dInit(Grid *g) {
  
  int choice;
  
  printf("Do you want 2D snapshots of the 3D grid? (1=yes, 0=no) ");
  scanf("%d", &choice);
  if (choice == 0) {
    temporalStride = -1;
    return;
  }
  
  printf("Duration of simulation is %d steps.\n", MaxTime);
  printf("Enter start time and temporal stride: ");
  scanf(" %d %d", &startTime, &temporalStride);
  printf("Enter the base name: ");
  scanf(" %s", basename);

  return;
}  /* end snapshot3dInit() */


void snapshot3d(Grid *g) {
  int mm, nn, pp;
  float temp;
  char filename[100];
  FILE *out;
  
  /* ensure temporal stride set to a reasonable value */
  if (temporalStride == -1) {
    return;
  } if (temporalStride < -1) {
    fprintf(stderr,
      "snapshot2d: snapshotInit2d must be called before snapshot.\n"
      "            Temporal stride must be set to positive value.\n");
    exit(-1);
  }

  /* get snapshot if temporal conditions met */
  if (Time >= startTime && 
      (Time - startTime) % temporalStride == 0) {

    /************ E slice ************/
    // printf("Hi1");
    sprintf(filename, "Data/%sE.%d.dat", basename, frameX++);
    out = fopen(filename, "w");
    
     /* write remaining data */
    pp = 1;
    for (nn = 0; nn < SizeY; nn++)
      for (mm = 0; mm < SizeX; mm++) {
	temp = (float)(sqrt(creal(Ex(mm, nn, pp))*creal(Ex(mm, nn, pp)) + cimag(Ex(mm, nn, pp)*cimag(Ex(mm, nn, pp))))); // store data as a float
	fprintf(out, "%f  \n", temp);
      }
    
    fclose(out);  // close file

    /**H slice**/
    //   printf("Hi2");
    sprintf(filename, "Data/%sH.%d.dat", basename, frameY++);
    out = fopen(filename, "w");
    
     /* write remaining data */
    pp = 1;
    for (nn = 0; nn < SizeY; nn++)
      for (mm = 0; mm < SizeX; mm++) {
	temp = (float)(sqrt(creal(Hz(mm, nn, pp))*creal(Hz(mm, nn, pp)) + cimag(Hz(mm, nn, pp)*cimag(Hz(mm, nn, pp))))); // store data as a float
	fprintf(out, "%f  \n", temp);
      }
    
    fclose(out);  // close file

    
    /************ write the constant-X slice ************/
  /*       sprintf(filename, "Data/%sZY.%d.dat", basename, frameY++); */
  /*   out = fopen(filename, "w"); */

    

  /*   /\* write remaining data *\/ */
  /*   mm = SizeX / 2; */
  /*   for (pp = SizeZ - 1; pp >= 0; pp--) */
  /*     for (nn = 0; nn < SizeX - 1; nn++) { */
  /* 	temp = (float)Ex(mm, nn, pp); // store data as a float */
  /* 	fprintf(out, "%f  \n", temp); // write the float */
  /*     } */

  /*   fclose(out);  // close file */
    
   }
 
  return;
}  /* end snapshot3d() */
