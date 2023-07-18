/************************************************************************
Cold_Collapse.c
Max Coy, 7/23
Simulates collisionless cold collapse code of uniform random sphere
 ***********************************************************************/


#define NMAX 2048 /*Max # of particles*/
#include "cc.h"




int main(void){
   
   /*Seeding Pseudo-rng*/
   // srand(time(NULL));

   int n = 1024; /*# of particles*/
   double m[NMAX], x[NMAX][3], v[NMAX][3], a[NMAX][3]; /*Particle Data*/
   
   double e, e_init, r_v, eps2, dt, T;
   
   double dt_exp;
   printf("Enter dt step size: 2^");
   int check_input = scanf("%lf",&dt_exp);
   if(check_input < 1)
   {
      printf("Error during data input\n");
      return -1;
   }
   dt = pow(2.0,dt_exp); /*Time step*/
   eps2 = pow(2.0,-10); /*Potential Softening Parameter^2*/

   T = 1.; /*Total time to be integrated*/
   

   mkrsphere(n, m, x, v, &r_v, eps2);

   /*Finding Initial Energy to be used for stability analysis*/
   calc_energy(&e_init, n, m, x, v, eps2, &r_v);/*Technically could do this in m_s_df*/
   
   
   /*Setting Initial Accelerations of particles*/
   calc_force(n, m, x, a, eps2);

   /*Setting up output file to save data*/
   // FILE *stream;
   // stream = fopen("Cloud_Collapse_dt9_1.csv","w");

   /*Stepping through Simulation*/
   double index = 0.;
   do{
      /*Integration*/
      leap_frog(n, m, x, v, a, dt, eps2);

      /*Outputing Current State to file*/
      // for(int i = 0; i < n; i++)
      // {
      //    /*For convenience, only saving position variables*/
      //    fprintf(stream, "%1.12f %1.12f %1.12f,", x[i][0], x[i][1], x[i][2]);
      // }
      calc_energy(&e,n, m, x, v, eps2, &r_v);
      
      // fprintf(stream,"%1.12f %1.12f\n", e, r_v); /*Adding energy / virial ratio of step*/

      /*Advancing time*/
      index += 1;
   }while(index*dt < T);
   
   calc_energy(&e,n, m, x, v, eps2, &r_v);
   
   /*Adding some final calculations / identifying variables to output file */
   // fprintf(stream,"%1.12f, %1.12f, %1.12f, %1.12f, %1.12f, %1.12f", e_init, e, dt, index, eps2, T);
   // fclose(stream); /*Closing output file*/
   printf("Integrated time: %1.2f\n",T);
   printf("Error: %1.6e\n",fabs( (e_init - e) / e_init));
   return 0;
}

