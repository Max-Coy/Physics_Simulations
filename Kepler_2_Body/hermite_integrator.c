/*4th order Hermite Integrator for N body problem
main() function simulates simple 2 body kepler problem 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NMAX 256 /*Set small because this particular code is only doing small preliminary testing*/
#define sqr(x) ((x) * (x))

#define G 1
#define PI 3.14159


/*Purpose: Calculate Current Energy of the system
Inputs: energy, number of particles, mass array, position array,
velocity array, epsilon^2, virial ratio
*/
void calc_energy(double *, int, double[], double[][3], double[][3], double, double *);


/*Individual timestep scheme hermite integrator, returns the minimum particle time after integration
Inputs: number of particles, mass arry, position array, velocity array, acceleration array,
jerk array, snap array, crackle array, timestep array, time array, softening parameter, 
accuracy parameter eta, minimum time (for output), fix_dt (inputting 0 means variable timestep is used)
*/
void its_hermite_integrator(int, double [], double [][3], double [][3], double [][3], 
                            double [][3], double [][3], double [][3], double [], double [], 
                            double, double, double *, int);

/*Function to be called prior to first call of its_hermite_integrator, sets initial t, dt, acceleration, and jerk values
Inputs: number of particles, mass arry, position array, velocity array, acceleration array,
jerk array, timestep array, time array, softening parameter, initial accuracy parameter eta_s, 
minimum allowed timestep, maximum allowed timestep, fixed_timestep (entering nonzero value will
set all dt values to this instead of doing calculation)
Note: Integrator does not impose restrictions on allowed timesteps so the timestep range is only for initialization
*/
void initialize_parameters(int, double [], double [][3], double [][3], double [][3], 
                            double [][3], double [][3], double [][3], double [], double [], 
                            double, double, int, int, double);

/*Finds specified time and outputs particle data at that point (right now just position and velocity)
Inputs: position output array, velocity output array, specified time, number of particles, position array, 
velocity array,  acceleration array, jerk array, snap array, crackle array, time array
WARNING: Onus is on user to make sure output_time is reasonable,
*/
void its_output(double [][3], double [][3], double, int, double [][3], double [][3], 
                double [][3], double [][3], double [][3], double [][3], double []);

/*Going to write a simple 2-body kepler problem in order to run preliminary analysis on stablity of 
integrator as well as ensure that things are working properly, as it won't be the meat of what's running
i'm not going to bother with small optimizations like avoiding divisions and such*/
int main(void)
{
    /*Initializing parameters*/
    double ms = 1e7, mp = 1e4; /*Mass of sun and planet respectively*/
    double sep = 20; /*initial separation of objects*/
    double rel_COM = sep*mp / (ms + mp); /*Finding relative center of mass bewteen objects (assuming sun is at [0,0] )*/
    /*Setting initial x-y coordinates of planet and sun so that center of mass is at [0,0]*/
    double rs[3] = {-rel_COM, 0, 0}; /*including z-coordinate so that dimensions match up with what code expects*/
    double rp[3] = {sep - rel_COM, 0, 0}; 

    /*Finding sufficient angular velocity for objects to ensure orbiting motion*/
    double eps2 = pow(2,-10); /*Potential softening parameter*/
    double fg = G * ms * mp / (sqr(sep) + eps2); /*Force of gravity acting between objects*/
    double mu = (ms*mp) / (ms + mp); /*reduced mass*/
    double omega = sqrt(fg / (mu * rp[0])); /*angular velocity (calculated around point of rotation and so rp instead of sep is used)*/
    
    /*Converting angular velocity to initial velocities*/
    double vs[3] = {-rs[1] * omega, rs[0] * omega, 0}; /*velocity of sun */
    double vp[3] = {-rp[1] * omega, rp[0] * omega, 0}; /*velocity of planet */

    /*Determining period of orbit*/
    double period = sqrt(4 * sqr(PI) * pow((sep * 0.5),3) / (G * (ms+mp) )  );
    
    /*Converting parameters into format the integrator can use*/
    double x[NMAX][3], v[NMAX][3], a[NMAX][3], j[NMAX][3], s[NMAX][3], c[NMAX][3]; /*Position and its derivatives*/
    double output_x[NMAX][3], output_v[NMAX][3]; /*Temp arrays to be used for outputs*/
    double eta = 0.005, eta_s = 0.001; /*Control accuracy / timestep size */
    int n = 2; /*number of particles*/
    double m[2] = {ms, mp}; /*particle mass array*/
    for(int i = 0; i < 3; i++) /*Setting initial position and velocity to previously calculated values*/
    {
        x[0][i] = rs[i];
        x[1][i] = rp[i];
        v[0][i] = vs[i];
        v[1][i] = vp[i];
    }
    double t[NMAX], dt[NMAX]; /*arrays to hold particle times and timesteps respectively*/
    double energy, energy0, r_v, t_min;
    double Time = 30 * period * 2.84423; /*Total time to be integrated (period seems to be off by a factor so I'm empirically adjusting it)*/
    t_min = 0.0;
    double index = 0.0;
    double resolution = pow(2,-12); /*we might have timesteps larger than this, will have to check*/

    double ts = pow(2,-16);

    calc_energy(&energy0, n, m, x, v, eps2, &r_v);

    initialize_parameters(n, m, x, v, a, j, s, c, dt, t, eps2, eta_s, 4, 16, 0.0);

    /*Outputting initial state*/
    FILE *stream;
    stream = fopen("hermite.csv", "w");

    for(int i = 0; i < n; i++)
    {
        fprintf(stream, "%1.12f %1.12f %1.12f ", x[i][0], x[i][1], x[i][2]);
        fprintf(stream, "%1.12f %1.12f %1.12f,", v[i][0], v[i][1], v[i][2]);
    }
    fprintf(stream, "%1.12f %1.12f\n",energy0,r_v);

    do
    {
        its_hermite_integrator(n, m, x, v, a, j, s, c, dt, t, eps2, eta, &t_min, 0);
        if(t_min > index * resolution)
        {
            its_output(output_x, output_v, index*resolution, n, x, v, a, j, s, c, t);
            calc_energy(&energy, n, m, output_x, output_v, eps2, &r_v);
            for(int i = 0; i < n; i++)
            {
                fprintf(stream, "%1.12f %1.12f %1.12f ", output_x[i][0], output_x[i][1], output_x[i][2]);
                fprintf(stream, "%1.12f %1.12f %1.12f,", output_v[i][0], output_v[i][1], output_v[i][2]);
            }
            fprintf(stream, "%1.16f %1.12f\n",energy,r_v);
            index++;
        }
    }while(index * resolution < Time);


    return 0;
}

/*Individual timestep scheme hermite integrator
Inputs: number of particles, mass arry, position array, velocity array, acceleration array,
jerk array, snap array, crackle array, timestep array, time array, softening parameter, 
accuracy parameter eta, minimum time (for output)
*/
void its_hermite_integrator(int n, double m[], double x[][3], double v[][3], double a[][3],
                        double j[][3], double s[][3], double c[][3], double dt[], double t[], 
                        double eps2, double eta, double *t_min, int fix_dt)
{
    /*Looping through particles to find particles with minimum t_i + dt_i*/
    double global_time; /*Variable to hold minimum time*/
    int indii[NMAX]; /*Array to hold indexes of particles*/
    int indexes = 1; /*Counter to keep track used length of array*/

    indii[0] = 0; /*Initializing minimum time as values for particle at index 0*/
    global_time = dt[0] + t[0];
    for(int i = 1; i < n; i++)
    {
        if((dt[i] + t[i]) < global_time)/*New global minimum found*/
        {
            global_time = dt[i] + t[i]; /*Setting new global minimum value*/
            indii[0] = i; /*Restarting array of indexes*/
            indexes = 1; /*Note that values past indexes will hold onto old/erraneous information*/
        }
        else if((dt[i] + t[i]) == global_time)/*Particle with degenerate time found*/
        {
            indii[indexes] = i; /*Adding to array of indexes*/
            indexes++; /*Increasing used length of array*/
        }
    }

    /*Predicting positions / velocities of all particles at global time*/
    double x_pred [NMAX][3], v_pred [NMAX][3]; /*Arrays to hold temporary data*/
    double t_factor, t_factor2, t_factor3; /*Temporary variables to save computational time*/
    double inv3 = 1.0/3.0;
    for(int i = 0; i < n; i++)
    {
        /*Predefining prefactors to reduce computations*/
        t_factor = global_time - t[i]; /* delta t */
        t_factor2 = sqr(t_factor) * 0.5; /* delta t^2 / 2 */
        t_factor3 = t_factor * t_factor2 * inv3; /* delta t^3 / 6 */

        for(int k = 0; k < 3; k++) /* j is already being used for jerk array */
        {
            x_pred[i][k] = t_factor3*j[i][k] + t_factor2*a[i][k] + t_factor*v[i][k] + x[i][k];
            v_pred[i][k] = t_factor2*j[i][k] + t_factor*a[i][k] + v[i][k];
        }
    }

    
     /*Setting array values for acceleration and jerk for all degenerate particles to 0 and storing them in temporary variables*/
    double a0[NMAX][3], j0[NMAX][3];
    int index;
    for(int i = 0; i < indexes; i++)
    {
        index = indii[i]; /*Using temporary variable to avoid excessive calls to potentially large array*/
        for(int k = 0; k < 3; k++)
        {
            a0[index][k] = a[index][k];
            j0[index][k] = j[index][k];

            a[index][k] = 0.0;
            j[index][k] = 0.0;
        }
    }
    
    /*Finding new acceleration and jerk for degenerate particles*/
    double dr[3], dv[3];
    double vDotR, R, rInv, rInv3, rInv5;

    for(int i = 0; i < indexes; i++)
    {
        index = indii[i]; /*Using temporary variable to avoid excessive calls to potentially large array*/
        for(int k = 0; k < n; k++)
        {
            if(k == index)
            {
                continue; /*Skipping this iteration as it's looking at self-interaction*/
            }
            vDotR = 0.0;
            /*Finding displacements*/
            for(int l = 0; l < 3; l++)
            {
                dr[l] = x_pred[k][l] - x_pred[index][l];
                dv[l] = v_pred[k][l] - v_pred[index][l];
                vDotR += dr[l] * dv[l];
            }
            /*Predefining prefactors for computations*/
            vDotR *= 3.0;
            R = sqrt(sqr(dr[0]) + sqr(dr[1]) + sqr(dr[2]) + eps2);
            rInv = 1 / R;
            rInv3 = sqr(rInv) * rInv;
            rInv5 = rInv3 * sqr(rInv);
            for(int l = 0; l < 3; l++)
            {
                a[index][l] += m[k] * dr[l] * rInv3;
                j[index][l] += m[k] * dv[l] * rInv3 + vDotR * dr[k] * rInv5;
            }
            /*Technically some calculations are being done (but not counted!) twice, specifically between any particles
            noted on the array indii, couldn't think of a clever way to avoid this cleanly right now but could but an
            array for some speed improvement*/
        }
    }


    /*Calculating corrections to predictions of degenerate particles and updating time*/
    double DT, dtI, dtI2, dtI3;
    double dt3, dt4, dt5;
    double inv_6 = 1.0 / 6.0;
    double inv_24 = inv_6 * 0.25;
    double inv_120 = inv_24 * 0.2;
    for(int i = 0; i < indexes; i++)
    {
        /*Calculating snap and crackle to use in corrections to predicted location of particle_index*/
        index = indii[i]; /*Using temporary variable to avoid excessive calls to potentially large array*/
        DT = dt[index];/*Not given that all dt's are the same so need to calc prefactors for each particle*/
        dt3 = sqr(DT) * DT;
        dt4 = DT * dt3;
        dt5 = DT * dt4;
        dtI = 1 / DT;
        dtI2 = sqr(dtI);
        dtI3 = dtI * dtI2;
        for(int k = 0; k < 3; k++)
        {
            /*Computing higher order derivatives*/
            s[index][k] = (-6 * (a0[index][k] - a[index][k]) - DT*(4 * j0[index][k] + 2*j[index][k]) ) * dtI2; /*Saving values to be used for timestep adjustment*/
            c[index][k] = (12 * (a0[index][k] - a[index][k]) + 6*DT*(j0[index][k] + j[index][k]) ) * dtI3; /*Saving values to be used for timestep adjustment*/

            /*Adjusting predictions and casting values to actual array*/
            x[index][k] = x_pred[index][k] + dt4 * inv_24 * s[index][k] + dt5 * inv_120 * c[index][k];
            v[index][k] = v_pred[index][k] + dt3 * inv_6  * s[index][k] + dt4 * inv_24  * c[index][k];
        }

        /*Updating particle time*/
        t[index] += dt[index];
    }

    /*Updating timesteps of degenerate particles*/
    if(fix_dt == 0) /*Input specifies that we want code to use variable timestep*/
    {
        double accel, jerk, jerk2, snap, snap2, crackle, dtNew;
        for(int i = 0; i < indexes; i++)
        {
            index = indii[i];

            /*First finding updated value for snap*/
            for(int k = 0; k < 3; k++)
            {
                s[index][k] += c[index][k] * dt[index]; /*no need to adjust crackle value b/c using 3rd order interpolation*/
            }

            /*Finding magnitudes of acceleration, jerk, snap, and crackle*/
            accel = sqrt( sqr(a[index][0]) + sqr(a[index][1]) + sqr(a[index][2]) );
            jerk2 = sqr(j[index][0]) + sqr(j[index][1]) + sqr(j[index][2]);
            jerk = sqrt(jerk2);
            snap2 = sqr(s[index][0]) + sqr(s[index][1]) + sqr(s[index][2]);
            snap = sqrt(snap2);
            crackle = sqrt( sqr(c[index][0]) + sqr(c[index][1]) + sqr(c[index][2]) ); 

            /*Calculating candidate timestep*/
            dtNew = sqrt(eta * (accel * snap + jerk2) / (jerk * crackle + snap2));

            if(dtNew < dt[index])
            {
                dt[index] *= 0.5;
            }
            else if(dtNew > 2 *dt[index]) /*dt[index] remains unchanged if dt[index] < dtNew <= 2*dt[index]*/
            {
                dt[index] *= 2;
            }
        }
    }

    /*Getting minimum particle time*/
    *t_min = t[0];
    for(int i = 1; i < n; i++)
    {
        // /*Converting following logic into branchless code
        if (t[i] < *t_min )
        {
            *t_min = t[i];
        }
        // */
    //    *t_min = (*t_min * (double) (*t_min <= t[i])) + (t[i] * (double) (*t_min > t[i]));
    }

    return;
}

/*Function to be called prior to first call of its_hermite_integrator, sets initial t, dt, acceleration,  
jerk, snap and crackle values
Inputs: number of particles, mass arry, position array, velocity array, acceleration array,
jerk array, timestep array, time array, softening parameter, initial accuracy parameter eta_s, 
*/
void initialize_parameters(int n, double m[], double x[][3], double v[][3], double a[][3],
                        double j[][3], double s[][3], double c[][3], double dt[], double t[], 
                        double eps2, double eta_s, int min_pow, int max_pow, double dt_fix)
{   
    double dr[3], dv[3];
    double vDotR, R, rInv, rInv3, rInv5, weight;

    /*Initializing t, dt, a, j, s, and c values to 0*/
    for(int i = 0; i < n; i++)
    {
        for(int k = 0; k < 3; k++)
        {
            a[i][k] = 0.0;
            j[i][k] = 0.0;
            s[i][k] = 0.0;
            c[i][k] = 0.0;
        }
        t[i] = 0.0;
        dt[i] = 0.0;
    }

    /*Setting initial a and j values*/
    for(int i = 0; i < n-1; i++)
    {
        for(int k = i + 1; k < n; k++)
        {
            vDotR = 0.0;
            /*Finding displacements*/
            for(int l = 0; l < 3; l++)
            {
                dr[l] = x[k][l] - x[i][l];
                dv[l] = v[k][l] - v[i][l];
                vDotR += dr[l] * dv[l];
            }
            /*Predefining prefactors for computations*/
            vDotR *= 3.0;
            R = sqrt(sqr(dr[0]) + sqr(dr[1]) + sqr(dr[2]) + eps2);
            rInv = 1 / R;
            rInv3 = sqr(rInv) * rInv;
            rInv5 = rInv3 * sqr(rInv);
            for(int l = 0; l < 3; l++)
            {
                weight = dr[l] * rInv3;
                a[i][l] += m[k] * weight;
                a[k][l] -= m[i] * weight;

                weight = dv[l] * rInv3 + vDotR * dr[l] * rInv5;
                j[i][l] += m[k] * weight;
                j[k][l] -= m[i] * weight;
            }
        }
    }
    
    /*Creating array to hold powers of 2 in specified range*/
    if((max_pow - min_pow + 1) > 32) /*Ensuring firs that we don't get an out of bound index*/
    {
        printf("Range of powers specified is too large, unable to set dt values");
        return;
    }
    double powers[32]; /*Can't set variable length so hard coding length of 32 which should be overkill*/
    int index = 0;
    for(int i = min_pow; i < (max_pow+1); i++)
    {
        powers[index] = pow(2.0, - (double) i);
        index++;
    }

    /*Setting initial dt values*/
    if(dt_fix == 0.0)
    {
        double unbound_dt;
        for(int i = 0; i < n; i++) //eta_s = 0.01 is supposedly a good value
        {
            unbound_dt = eta_s * sqrt(sqr(a[i][0]) + sqr(a[i][1]) + sqr(a[i][2]))
                        / sqrt(sqr(j[i][0]) + sqr(j[i][1]) + sqr(j[i][2]));
            
            /*Flooring value to next power of 2*/
            index = max_pow - min_pow;
            if(unbound_dt <= powers[max_pow - min_pow]) /*Checking to see if value is below allowed range*/
            {
                dt[i] = powers[max_pow - min_pow]; /*Setting to min value*/
                index = 0;/*Prevents us from entering while-loop*/
            }else if(unbound_dt >= powers[0]) /*Checking to see if above allowed range*/
            {
                dt[i] = powers[0]; /*Setting to max value*/
                index = 0;/*Prevents us from entering while-loop*/
            }
            while(index > 0) /*Looping through list from smallest to highest values*/
            {
                index--; /*already checked smallest value so first we index down*/
                if(unbound_dt == powers[index])
                {
                    dt[i] = powers[index];
                    break;
                }else if(unbound_dt < powers[index])
                {
                    dt[i] = powers[index+1];
                    break;
                }
            }
        }
    }else{/*dt value was specifed for initialization*/
        for(int i = 0; i < n; i++)
        {
            dt[i] = dt_fix;
        }
    }
    return;
}

/*Finds specified time and outputs particle data at that point (right now just position)
Inputs: output array, specified time, number of particles, position array, velocity array, 
acceleration array, jerk array, snap array, crackle array, time array
WARNING: Onus is on user to make sure output_time is reasonable,
*/
void its_output(double x_pred[][3], double v_pred[][3], double output_time, int n, double x[][3],  
                double v[][3], double a[][3],  double j[][3], double s[][3], double c[][3], double t[])
{    
    /*Temporary variables to save computational time*/
    double dt1, dt2, dt3, dt4, dt5; 
    double inv_3 = 1.0 / 3.0;

    /*Predicting positions / velocities of all particles at global time*/
    for(int i = 0; i < n; i++)
    {
        /*Predefining prefactors to reduce computations*/
        dt1 = output_time - t[i]; /* delta t */
        dt2 = dt1 * dt1 * 0.5; /* delta t^2 / 2 */
        dt3 = dt1 * dt2 * inv_3; /* delta t^3 / 6 */
        dt4 = dt1 * dt3 * 0.25; /* delta t^4 / 24 */
        dt5 = dt1 * dt4 * 0.2; /* delta t^5 / 120 */
        for(int k = 0; k < 3; k++) /* j is already being used for jerk array */
        {
            x_pred[i][k] = x[i][k] + v[i][k] * dt1 + a[i][k] * dt2 + j[i][k] * dt3
                         + s[i][k] * dt4 + c[i][k] * dt5;
            v_pred[i][k] = v[i][k] + a[i][k] * dt1 + j[i][k] * dt2 + s[i][k] * dt3
                         + c[i][k] * dt4;
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
