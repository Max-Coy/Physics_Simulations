/************************************************************************
leap_frog.c
Max Coy, 7/23
Simple leap_frog integration scheme for N-body gravitational problem
 
Functions:
   calc_force()      Int      calculates inter-particle forces
   leap_frog()       Ext      leap_frog integrator
   calc_energy()     Ext      calculates system energy

External Variables:
   double   a[][]    I/O      accelerations of particles
   double   dt       I        timestep of simulation
   double   eps2     I        potential softening parameter
   double   m[]      I        mass of particles
   double   v[][]    I/O      velocity of particles
   double   x[][]    I/O      posiiton of particles

Requires functions from "mfunc.c" to run
***********************************************************************/


/*Inter-particle force calculator to set accelerations*/
void calc_force(int n, double m[], double x[][3], double a[][3], double eps2)
{
   /*Clearing / Initializing Acceleration Data*/
   for(int i = 0; i < n; i++)
   {
      for(int j = 0; j < 3; j++)
      {
         a[i][j] = 0.0;
      }
   }

   double r[3]; /*array to hold relative cordinate displacements*/
   double R2, R3_inv; /*variable to hold Total displacement*/

   /*Calculating Accelerations*/
   for(int i = 0; i < n-1; i++)
   {
      for(int j = i+1; j < n; j++)
      {
         /*Finding Displacements*/
         for(int k = 0; k < 3; k++)
         {
            r[k] = x[j][k] - x[i][k]; 
         }


         R2 = mag3_sqr(r) + eps2;
         R3_inv = sqrt(1 / (R2*R2*R2)); /*Inverse cube of R*/

         /*Adjusting Accelerations*/
         for(int k = 0; k < 3; k++)
         {
            a[i][k] +=  m[j] * r[k] * R3_inv;
            a[j][k] += -m[i] * r[k] * R3_inv;
         }
      }
   }
   return;
}

/*Leapfrog Numerical Integrator*/
void leap_frog(int n, double m[], double x[][3], double v[][3], double a[][3],
               double dt, double eps2)
{
   
   double dt_half = dt * 0.5; /*saves computation time*/

   /*Adjusting Positions*/
   for(int i = 0; i < n; i++)
   {
      for(int j = 0; j < 3; j++)
      {
         v[i][j] += a[i][j]*dt_half; /* Going 1/2 timestep ahead for v*/
         x[i][j] += v[i][j]*dt;
      }
   }

   /*Adjusting Acclerations based on new positions*/
   calc_force(n, m, x, a, eps2);


   /*Taking 2nd half-step for velocities*/
   for(int i = 0; i < n; i++)
   {
      for(int j = 0; j < 3; j++)
      {
         v[i][j] += a[i][j] * dt_half;
      }
   }

   return;
}


/*Calculates current total energy of system*/
void calc_energy(double *e, int n, double m[], double x[][3], double v[][3], 
                 double eps2, double *r_v)
{
   double K = 0.0, W = 0.0; /*Kinetic and Potential Energies*/
   double r2;
   
   /*Calculating Kinetic Energy*/
   for(int i = 0; i < n; i++)
   {
      K += m[i] * mag3_sqr(v[i]); 
   }
   K *= 0.5; /*Adding missing 1/2 factor */

   /*Calculating Potential Energy*/
   for(int i = 0; i < n-1; i++)
   {
      for(int j = i+1; j < n; j++)
      {
         r2 = dist3_sqr(x[j], x[i]);
         W -= m[i]*m[j] / sqrt(r2 + eps2);
      }
   }

   /*Calculating Virial Ratio*/
   *r_v = K/ fabs(W);
  
   *e = (K + W);

   return;
}
