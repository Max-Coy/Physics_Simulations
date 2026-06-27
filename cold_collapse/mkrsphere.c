/******************************************************************************
mkrsphere.c
Max Coy, 7/23
Generates a uniform random spherical distribution in cartesian coordinates

Functions:
    guassian()          Int     Guassian distribution 
    mkrsphere()         Ext     Generates uniform random sphere

External Variables:
    double  eps2        I       potential softening parameter
    double  m[]         O       mass of particles
    int     M           I       total mass of system
    int     n           I       total number of particles
    int     NMAX        SC      maximum allowed number of particles
    double  r_v         O       virial ratio
    double  v[][]       O       velocity of particles
    double  x[][]       O       position of particles

Requires functions from "mfunc.c" to run
******************************************************************************/


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

/*Sets up initial conditions for simulation*/
void mkrsphere(int n, double m[], double x[][3], double v[][3],
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
         r2 = dist3_sqr(x[j], x[i]);
         W -= m[i]*m[j] / sqrt(r2 + eps2);
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


