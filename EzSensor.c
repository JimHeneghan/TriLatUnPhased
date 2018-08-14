#include <stdio.h>
#include <stdlib.h>
#include "fdtd-macro.h"
#include <complex.h>
#include <math.h>
/*code to record a sensor of the data*/
void sensorInit(Grid *g, double Kx, double Ky, double Kz, int j){
	char filename[100];
	FILE *out;
	sprintf(filename, "PBGData/Aug9Phased/EmptyLat%d.txt",j);
        out = fopen(filename, "a");
	fprintf(out, "%g \t %g \t %g \n", Kx, Ky, Kz);
	fclose(out);

return;
}



void Transmission(Grid *g, int j,  double time){
        char filename[100];
        double Time1;
        FILE *out;
        sprintf(filename, "PBGData/Zone6/EmptyLat%d.txt", j);
        out = fopen(filename, "a");
        Time1 = time*5.4683e-13;

        /* print the time stamp and the Ex field right before the QWS*/
        fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n", Time1, creal(Hz(80, 34, 1)), cimag(Hz(80, 34, 1)), creal(Hz(33, 43, 1)), cimag(Hz(33, 43, 1)), creal(Hz(83, 8, 1)), cimag(Hz(83, 8, 1)), creal(Hz(81, 44, 1)), cimag(Hz(81, 44, 1)), creal(Hz(36, 44, 1)), cimag(Hz(36, 44, 1)));
        fclose(out);

        char filename1[100];
        sprintf(filename1, "FishData/Ex%d.txt", j);
        out = fopen(filename1, "a");
        Time1 = time*5.4683e-13;

        /* print the time stamp and the Ex field right before the QWS*/
        fprintf(out, "%g \t %g \t %g \n", Time1, creal(Ex(80, 34, 1)), cimag(Ex(80, 34, 1)));
        fclose(out);

        char filename2[100];
        sprintf(filename2, "FishData/Ey%d.txt", j);
        out = fopen(filename2, "a");
        Time1 = time*5.4683e-13;

        /* print the time stamp and the Ex field right before the QWS*/
        fprintf(out, "%g \t %g \t %g \n", Time1, creal(Ey(80, 34, 1)), cimag(Ey(80, 34, 1)));
        fclose(out);

        char filename3[100];
        sprintf(filename3, "FishData/Ez%d.txt", j);
        out = fopen(filename3, "a");
        Time1 = time*5.4683e-13;

        /* print the time stamp and the Ex field right before the QWS*/
        fprintf(out, "%g \t %g \t %g \n", Time1, creal(Ez(80, 34, 1)), cimag(Ez(80, 34, 1)));
        fclose(out);

        char filename6[100];
        sprintf(filename6, "FishData/Hx%d.txt", j);
        out = fopen(filename6, "a");
        Time1 = time*5.4683e-13;

        /* print the time stamp and the Ex field right before the QWS*/
        fprintf(out, "%g \t %g \t %g \n", Time1, creal(Hx(80, 34, 1)), cimag(Hx(80, 34, 1)));
        fclose(out);

        char filename4[100];
        sprintf(filename4, "FishData/Hy%d.txt", j);
        out = fopen(filename4, "a");
        Time1 = time*5.4683e-13;

        /* print the time stamp and the Ex field right before the QWS*/
        fprintf(out, "%g \t %g \t %g \n", Time1, creal(Hy(80, 34, 1)), cimag(Hy(80, 34, 1)));
        fclose(out);

        char filename5[100];
        sprintf(filename5, "FishData/Hz%d.txt", j);
        out = fopen(filename5, "a");
        Time1 = time*5.4683e-13;

        /* print the time stamp and the Ex field right before the QWS*/
        fprintf(out, "%g \t %g \t %g \n", Time1, creal(Hz(80, 34, 1)), cimag(Hz(80, 34, 1)));
        fclose(out);
}



