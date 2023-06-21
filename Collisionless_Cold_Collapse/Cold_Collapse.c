/******************************************************************************
 Collisionless Cold Collapse
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define NMAX 2048 /*Max # of particles*/
//#define std 1 /*STD used w/ guassian function*/
#define SQR(x) ((x)*(x))
#define drand48() ((double)rand() / RAND_MAX )

#define G 1 /* Gravitational Constant*/
#define M 1 /* Total Mass of system*/


/*
Purpose: Set Initial Conditions of particles
Inputs: Number of particles, mass array, position array, 
velocity array, virial ratio, epsilon^2
*/
void make_spherical_df(int, double[], double [][3], double [][3],
                       double *, double);

/* Guass with mean = 0.0 and dispersion = 1.0 by Box-Muller*/
double guassian(void);

/*
Purpose: Calculates inter-particle forces to modify acceleration array
Inputs: number of particles, mass array, position array,
acceleration array, epsilon^2
*/
void calc_force(int, double[], double[][3], double[][3], double);

/*
Purpose: Numerical Integrator
Inputs: number of particles, mass array, position array,
velocity array, acceleration array, timestep, epsilon
*/
void leap_frog(int, double[], double[][3], double[][3], double[][3],
               double, double);

/*
Purpose: Calculate Current Energy of the system
Inputs: energy, number of particles, mass array, position array,
velocity array, epsilon^2, virial ratio
*/
void calc_energy(double *, int, double[], double[][3], double[][3], double, double *);

int main(void){
   
   /*Seeding Pseudo-rng*/
   srand(time(NULL));



   int n = 1024; /*# of particles*/
   double m[NMAX], x[NMAX][3], v[NMAX][3], a[NMAX][3]; /*Particle Data*/
   
   double e, e_init, r_v, eps2, dt, T;
   
   

   dt = pow(2.0,-5); /*Time step*/
   eps2 = pow(2.0,-10); /*Potential Softening Parameter^2*/

   T = 10.; /*Total time to be integrated*/
   

   make_spherical_df(n, m, x, v, &r_v, eps2);

   /*Finding Initial Energy to be used for stability analysis*/
   calc_energy(&e_init, n, m, x, v, eps2, &r_v);/*Technically could do this in m_s_df*/
   
   /*Setting Initial Accelerations of particles*/
   calc_force(n, m, x, a, eps2);

   /*Stepping through Simulation*/
   double index = 0.;
   do{
      /*Integration*/
      leap_frog(n, m, x, v, a, dt, eps2);

      calc_energy(&e,n, m, x, v, eps2, &r_v);
      
      /*Advancing time*/
      index += 1;
   }while(index*dt < T);
   
   calc_energy(&e,n, m, x, v, eps2, &r_v);
   
   return 0;
}


/*Numerical Integrator*/
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
   double R, R3I; /*variable to hold Total displacement*/

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

         R = sqrt(SQR(r[0]) + SQR(r[1]) + SQR(r[2]) + eps2); 
         R3I = 1 / (R*R*R); /*Inverse cube of R*/

         /*Adjusting Accelerations*/
         for(int k = 0; k < 3; k++)
         {
            a[i][k] +=  m[j] * r[k] * R3I;
            a[j][k] += -m[i] * r[k] * R3I;
         }
      }
   }
   return;
}


/*Sets up initial conditions for simulation*/
void make_spherical_df(int n, double m[], double x[][3], double v[][3],
                       double *r_v, double eps2)
{
   if (n > NMAX){ /*Ensuring that we are within limits of program*/
      printf("Error! number of particles exceeds allowed limit");
      return;
   }

   double mass = (double) M / (double) n; /*Individual Particle mass*/
   
   double x1, x2, x3, r2; /*To be used for setting initial positions*/
   double K = 0.0, W=0.0; /*Kinetic/Potential Energies, to be used for setting Virial Ratio*/
   double std; /*To be used for setting initial velocities*/

   for(int i = 0; i < n; i++)
   {
      /*Setting Mass*/
      m[i] = mass; 

      /*Setting Initial Position*/
      do{
         x1 = 2.0 * drand48() - 1.0; /*Generating random position*/
         x2 = 2.0 * drand48() - 1.0;
         x3 = 2.0 * drand48() - 1.0;
         r2 = x1*x1 + x2*x2 + x3*x3;
      }while (r2 > 1); /*Ensuring position lies within sphere*/

      x[i][0] = x1;
      x[i][1] = x2;
      x[i][2] = x3;
      
   }

   /*Calculating Potential Energy*/
   for(int i = 0; i < n-1; i++)
   {
      for(int j = i+1; j < n; j++)
      {
         r2 = SQR(x[j][0] - x[i][0]) + SQR(x[j][1] - x[i][1]) + SQR(x[j][2] - x[i][2]);
         W = W - m[i]*m[j] / sqrt(r2 + eps2);
      }
   }

   *r_v = 0.1; //setting initial value for r_v
   std = sqrt(2.0 * *r_v * fabs(W) / 3.0);


   /*Setting Initial Velocities*/
   for(int i = 0; i < n; i++)
   {
      for(int j = 0; j < 3; j++)
      {
         v[i][j] = std*guassian();
      }
   }



   return;
}

/*Calculates current total energy of system*/
void calc_energy(double *e, int n, double m[], double x[][3], double v[][3], 
                 double eps2, double *r_v)
{
   double K = 0.0, W = 0.0;
   double r2;
   
   /*Calculating Kinetic Energy*/
   for(int i = 0; i < n; i++)
   {
      K += m[i] * (SQR(v[i][0]) + SQR(v[i][1]) + SQR(v[i][2])); 
   }
   K *= 0.5; /*Adding missing 1/2 factor */

   /*Calculating Potential Energy*/
   for(int i = 0; i < n-1; i++)
   {
      for(int j = i+1; j < n; j++)
      {
         r2 = SQR(x[j][0] - x[i][0]) + SQR(x[j][1] - x[i][1]) + SQR(x[j][2] - x[i][2]);
         W -= m[i]*m[j] / sqrt(r2 + eps2);
      }
   }

   /*Calculating Virial Ratio*/
   *r_v = K/ fabs(W);
  
   *e = (K + W);

   return;
}


/* Guass with mean = 0.0 and dispersion = 1.0 by Box-Muller*/
double guassian(void)
{
   double x, y, r2, z;

   do{
      x = 2.0 * drand48() - 1.0;
      y = 2.0 * drand48() - 1.0;
      r2 = x*x + y*y;
   }while(r2 >= 1.0 || r2 == 0.0);

   z = sqrt(-2.0 * log(r2) / r2) * x;

   return z;
}
