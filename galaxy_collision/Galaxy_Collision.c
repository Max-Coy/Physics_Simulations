/****************************************
Galaxy_Collision.c
Max Coy, 7/23
Program to simulate galaxies colliding
****************************************/


#define NMAX 8192
#define M 1
#define G 1
#define DIM 3

#include "gc.h"


int main(void)
{
    double n1 = 1024, n2 = 512, n3 = 512; /*Respective particle counts for each cluster*/
    /*Ensuring that total particle count stays within allowed limits*/
    if (n1 + n2 > NMAX)
    {
        printf("Warning, total particle count exceeds allowed maximum, exiting code");
        return 0;
    }

    // srand(time(NULL));


    double m[NMAX], x[NMAX][3], v[NMAX][3], a[NMAX][3];
    double Energy, Energy_initial, r_v, eps2, dt, T;


    eps2 = pow(2,-10); /*Potential softening parameter*/
    dt = pow(2,-5); /*Timestep*/
    T = 0.5; /*Total time to be integrated*/

    /*Generating initial state*/
    double r_cut = 10; /*Cutoff radius for plummer sphere*/
    mkplummer(n1, m, x, v, r_cut);
    double xShift1[3] = {4, -2, 0};
    double vShift1[3] = {-0.25, 0.05, 0};
    shift_arr_values2(n1, 0, xShift1, x, DIM); /*Moving particles to right*/
    shift_arr_values2(n1,0,vShift1, v, DIM);/*Giving total cluster initial velocity*/
    shift_arr_indicies2(n1, 0, n2+n3, x, DIM); /*Making space in array to generate second cluster*/
    shift_arr_indicies2(n1, 0, n2+n3, v, DIM);
    shift_arr_indicies(n1, 0, n2+n3, m);

    mkplummer(n2, m, x, v, r_cut);
    double xShift2[3] = {-4, 2, 0};
    double vShift2[3] = {0.3, 0.125, 0};
    shift_arr_values2(n2, 0, xShift2, x, DIM); /*Shifting second cluser left*/
    shift_arr_values2(n2, 0, vShift2, v, DIM);
    shift_arr_indicies2(n2, 0, n3, x, DIM);
    shift_arr_indicies2(n2, 0, n3, v, DIM);
    shift_arr_indicies(n2, 0, n3, m);

    mkplummer(n3, m, x, v, r_cut);
    double xShift3[3] = {-5, -5, 0};
    double vShift3[3] = {-0.15, 0.3, 0};
    shift_arr_values2(n3, 0, xShift3, x, DIM); 
    shift_arr_values2(n3,0,vShift3, v, DIM);

    
    int ntot = n1 + n2 + n3;
    calc_energy(&Energy_initial, ntot, m, x, v, eps2, &r_v);
    calc_force(ntot, m, x, a, eps2);

    /*Setting up output file for data*/
    FILE *stream;
    stream = fopen("sanity.csv", "w");

    /*Stepping through sim*/
    double index = 0.;
    do{
      /*Integration*/
      leap_frog(ntot, m, x, v, a, dt, eps2);

      /*Outputing Current State to file*/
      for(int i = 0; i < ntot; i++)
      {
         /*For convenience, only saving position variables*/
         fprintf(stream, "%1.12f %1.12f %1.12f,", x[i][0], x[i][1], x[i][2]);
      }
      calc_energy(&Energy, ntot, m, x, v, eps2, &r_v);
      fprintf(stream,"%1.12f %1.12f\n",Energy,r_v);

      /*Advancing time*/
      index += 1;
    }while(index*dt < T);

    /*Outputting some final descriptive variables*/
    fprintf(stream,"%1.12f, %1.12f, %1.12f, %1.12f, %1.12f, %1.12f", Energy_initial, 
            Energy, dt, index, eps2, T);
    fclose(stream);

    return 0;
}
