/*Galaxy_Collision.c
Program to simulate galaxies colliding
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NMAX 8192
#define pi 3.14159
#define inv_pi 0.3183
#define M 1
#define G 1

#define sqr(x) ((x) * (x))
#define drand() ( (double)rand() / RAND_MAX )

/*Calculates the scale factor of the total mass included within the cutoff radius
Inputs: cutoff radius
*/
double calc_mass_cutoff(double);

/*Creates a radial distribution using the given mass cutoff factor
Inputs: number of particles, array for radial values, mass cutoff factor
*/
void radial_distribution(int, double [], double);

/*Generates starting positions of particles in cartesian coordinates
Inputs: number of particles, array for radial values, array for positions 
*/
void generate_particle_distribution(int, double [], double [][3]);

/*Generates initial velocities of particles, dependant on initial radius
Inputs: number of particles, velocity array, radial array
*/
void generate_velocity_distribution(int, double [][3], double [], double [][3]);

/* Initialize particles into plummer model distribution
Inputs: number of particles, mass array, position array, velocity array, softening parameter, virial radius
*/
void mkplummer(int, double [], double [][3], double [][3], double);

/*Purpose: Calculates inter-particle forces to modify acceleration array

Inputs: number of particles, mass array, position array,
acceleration array, epsilon^2
*/
void calc_force(int, double[], double[][3], double[][3], double);

/*Purpose: Numerical Integrator

Inputs: number of particles, mass array, position array,
velocity array, acceleration array, timestep, epsilon
*/
void leap_frog(int, double[], double[][3], double[][3], double[][3],
               double, double);

/*Purpose: Calculate Current Energy of the system
Inputs: energy, number of particles, mass array, position array,
velocity array, epsilon^2, virial ratio
*/
void calc_energy(double *, int, double[], double[][3], double[][3], double, double *);

/*Purpose: Adjust position of particle information within arrays
Inputs: Current offset, New offset, mass array, position array, velocity array
acceleration array, order of the last 3 inputs doesn't actually matter as it's 
an identical operation performed on all 3
*/
void shift_indicies(int, int, int, double [], double [][3], double [][3], double [][3]);

/*Purpose: Shift the specified values of particles by fixed amount
Inputs: number of particles, index offset, array containing desired offset,
particle position array
*/
void shift_values(int, int, double [3], double [][3]);

/*Purpose: Scale given array by given factor
Inputs: number of particles, current offset, scaling factor, array to be scaled
Note: not intended to be used for mass array
*/
void scale_values(int, int, double, double [][3]);


int main(void)
{
    double n1 = 1024, n2 = 512, n3 = 512; /*Respective particle counts for each cluster*/
    printf("a");
    /*Ensuring that total particle count stays within allowed limits*/
    if (n1 + n2 > NMAX)
    {
        printf("Warning, total particle count exceeds allowed maximum, exiting code");
        return 0;
    }

    srand(time(NULL));


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
    shift_values(n1, 0, xShift1, x); /*Moving particles to right*/
    shift_values(n1,0,vShift1, v);/*Giving total cluster initial velocity*/
    shift_indicies(n1, 0, n2+n3, m, x, v, a); /*Making space in array to generate second cluster*/
   
    mkplummer(n2, m, x, v, r_cut);
    double xShift2[3] = {-4, 2, 0};
    double vShift2[3] = {0.3, 0.125, 0};
    shift_values(n2, 0, xShift2, x); /*Shifting second cluser left*/
    shift_values(n2, 0, vShift2, v);
    shift_indicies(n2, 0, n3, m, x, v, a);

    mkplummer(n3, m, x, v, r_cut);
    double xShift3[3] = {-5, -5, 0};
    double vShift3[3] = {-0.15, 0.3, 0};
    shift_values(n3, 0, xShift3, x); 
    shift_values(n3,0,vShift3, v);

    
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



/*Stealing these functions from my cold collapse code with some tweaks*/
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

/*Force calculator to set accelerations*/
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

         R = sqrt(sqr(r[0]) + sqr(r[1]) + sqr(r[2]) + eps2); 
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

/*Calculates current total energy of system*/
void calc_energy(double *e, int n, double m[], double x[][3], double v[][3], 
                 double eps2, double *r_v)
{
   double K = 0.0, W = 0.0;
   double r2;
   
   /*Calculating Kinetic Energy*/
   for(int i = 0; i < n; i++)
   {
      K += m[i] * (sqr(v[i][0]) + sqr(v[i][1]) + sqr(v[i][2])); 
   }
   K *= 0.5; /*Adding missing 1/2 factor */

   /*Calculating Potential Energy*/
   for(int i = 0; i < n-1; i++)
   {
      for(int j = i+1; j < n; j++)
      {
         r2 = sqr(x[j][0] - x[i][0]) + sqr(x[j][1] - x[i][1]) + sqr(x[j][2] - x[i][2]);
         W -= m[i]*m[j] / sqrt(r2 + eps2);
      }
   }

   /*Calculating Virial Ratio*/
   *r_v = K/ fabs(W);
  
   *e = (K + W);

   return;
}


/*Upcoming is a few functions for shifting/scaling generated spheres*/
/*Purpose: Adjust position of particle information within arrays
Inputs: number of particles, Current offset, New offset, mass array, position array, 
velocity array acceleration array, order of the last 3 inputs doesn't actually matter 
as it's an identical operation performed on all 3
WARNING: This will write over data so make sure you don't botch the offsets
*/
void shift_indicies(int n, int current_offset, int new_offset, double m[], double x[][3], 
                    double v[][3], double a[][3])
    {
        /*Ensuring there is space in the array for new particles*/
        if(n + new_offset > NMAX)
        {
            printf("Warning, attempting to index variables out of array bounds, terminating call");
            return;
        }
        /*Adjusting information*/
        for(int i = 0; i < n; i++)
        {
            m[i+new_offset] = m[i+current_offset];
            for(int j = 0; j < 3; j++)
            {
                x[i+new_offset][j] = x[i+current_offset][j];
                v[i+new_offset][j] = v[i+current_offset][j];
                a[i+new_offset][j] = a[i+current_offset][j];
            } 
        }
        return;
    }

/*Purpose: Shift the specified values of particles by fixed amount
Inputs: number of particles, index offset, array containing desired offset,
particle value array
*/
void shift_values(int n, int offset, double shift[3], double arr[][3])
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            arr[i + offset][j] += shift[j];
        }
    }
    return;
}

/*Purpose: Scale given array by given factor
Inputs: number of particles, current offset, scaling factor, array to be scaled
Note: not intended to be used for mass array
*/
void scale_values(int n, int offset, double factor, double arr[][3])
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            arr[i + offset][j] *= factor;
        }
    }
    return;
}

/*Up Next is a series of functions used for generating the initial Plummer Sphere distribution*/
/*Calculates the scale factor of the total mass included within the cutoff radius
Inputs: cutoff radius
*/
double calc_mass_cutoff(double r_cut)
{
    /*In the case of large r, there is no appreciable mass loss*/
    if (r_cut > 99999.0)
    {
        return 1.0;
    }
    double scale_factor = 5.333 * inv_pi; /* 16 / 3 pi*/ 
    double r_scaled = scale_factor * r_cut;
    double num = r_scaled * r_scaled * r_scaled;
    double denom = pow(1.0 + r_scaled*r_scaled, 1.5);
    return num / denom; 

}

/*Creates a radial distribution using the given mass cutoff factor
Inputs: number of particles, array for radial values, mass cutoff factor
*/
void radial_distribution(int n, double r[], double mass_cut)
{
    double exp = -2.0/3.0; /*Exponent to be used later*/
    for(int i = 0; i < n; i++)
    {
        r[i] = sqrt(pow(drand()*mass_cut,exp) - 1); 
    }
    return;
}

/*Generates starting positions of particles in cartesian coordinates
Inputs: number of particles, array for radial values, array for positions 
*/
void generate_particle_distribution(int n, double r[], double x[][3])
{
    double theta, phi;
    for(int i = 0; i < n; i++)
    {
        theta = drand() * 2*pi;
        phi = acos(1 - drand()*2);
        x[i][0] = r[i] * sin(phi) * cos(theta);
        x[i][1] = r[i] * sin(phi) * sin(theta);
        x[i][2] = r[i] * cos(phi);
    }
    return;
}

/*Generates initial velocities of particles, dependant on initial radius
Inputs: number of particles, velocity array, radial array
*/
void generate_velocity_distribution(int n, double v[][3], double r[], double x[][3])
{
    /* Been doing some reading, if anything is going wrong with the code it'll probably be here, gonna follow solution
    outlined here: https://physics.stackexchange.com/questions/94845/velocity-distribution-in-plummers-models-and-others-mass-distributions
    while refrencing: https://github.com/amusecode/amuse/blob/main/src/amuse/ic/plummer.py 
    */
    double v_esc, weight, v_mag, y;
    double theta, phi, rho;
    double sqrt2 = sqrt(2);
    double v1[3], v2[3], vp[2];

    for(int i = 0; i < n; i++)
    {
        do
        {
            y = drand() * 0.1;
            weight = drand();
            weight = pow(1-sqr(weight),3.5) * sqr(weight);
        }while(y > weight); /*Half-understand / half don't understand what's happening in the loop lmao*/

        v_esc = sqrt2 * pow(1 + sqr(r[i]), -0.25); /*Escape velocity at given radius*/
        v_mag = v_esc*weight * 7; /*Total magnitude of particle velocity*/
        /*^ *7 is a fudge factor added to prevent collapse, not really sure what's going wrong with
        my velocity calc that is making me need it but alas...*/

        /*Pointing velocity in random direction*/
        theta = drand() * 2*pi;
        phi = acos(1 - drand()*2);

        v[i][0] = v_mag * sin(phi) * cos(theta);
        v[i][1] = v_mag * sin(phi) * sin(theta);
        v[i][2] = v_mag * cos(phi);

    }
    
    return;
}

/* Initialize particles into plummer model distribution
Inputs: number of particles, mass array, position array, velocity array, softening parameter, virial radius
*/
void mkplummer(int n, double m[], double x[][3], double v[][3], double r_cut)
{
    
    if(n > NMAX)
    {
        printf("Warning, particle count has exceeded allowed limit, exiting code");
        return;
    }

    double m_cut = calc_mass_cutoff(r_cut); /*Scale factor of mass cutoff*/
    double r[NMAX]; /*array holding radial values*/

    /*Setting all particles to have equal mass*/
    double mass = (double) M / (double) n;
    for(int i = 0; i < n; i++)
    {
        m[i] = mass;
    }

    /*Generating initial positions*/
    radial_distribution(n, r, m_cut);
    generate_particle_distribution(n, r, x);

    /*Generating initial velocities*/
    generate_velocity_distribution(n, v, r, x);

    return;
}
