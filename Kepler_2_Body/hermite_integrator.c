/****************************************************************************************************
hermite_integrator.c
Max Coy, 7/23
4th order (PEC)^n Hermite Integrator for N body gravitational problem 
Can be set to use variable timesets of integer powers of 2, or use a prespecified fixed timestep

Functions:
    its_hermite_integrator()    Ext         Individual Time Step hermite integrator
    initialize_parameteres()    Ext         initializes variable arrays for simulation
    its_output()                Ext         predicts particle X and V values for given time t
    calc_energy()               Ext         Calculates system energy

External Variables:
    double  a[][]               I/O         accelerations of particles
    double  c[][]               I/O         crackle of particles
    int     DIM                 SC          Total dimensions of space (just set it to 3)
    double  dt[]                I/O         timesteps of particles
    double  dt_fix              I           sets fixed initial timestep if nonzero
    double  e                   O           total energy of the system
    double  eps2                I           gravitational softening parameter
    double  eta                 I           timestep accuracy parameter
    double  eta_s               I           initial timestep accuracy parameter
    int     indexes             O           length of indii[]
    int     indii[]             O           indexes of particles adjusted by integrator
    int     fix_dt              I           disables individual setting of timesteps if true
    double  j[][]               I/O         jerk of particles
    double  m[]                 I/O         mass of particles
    double  max_pow             I           smallest initial timestep allowed (2^(-max_pow))
    double  min_pow             I           largest initial timestep allowed (2^(-min_pow))
    int     n                   I           number of particles in system
    int     NMAX                SC          max number of particles allowed in simulation
    int     order               I           order of the (PEC)^n scheme
    double  r_v                 O           virial ratio
    double  s[][]               I/O         snap of particles
    double  t[]                 I/O         time of particles
    double  t_min               O           minimum value in t[]
    double  v[][]               I/O         velocity of particles
    double  v_pred[][]          O           predicted velocity of particles
    double  x[][]               I/O         position of particles
    double  x_pred[][]          O           predicted position of particles

Requires functions from "mfunc.c" and "array_func.c" to run
Requires symbolic constants from "mconst.h" to run
WARNING: using -O2 seems to brick code so don't use it...
****************************************************************************************************/


/*Individual timestep scheme hermite integrator
Inputs: number of particles, order of scheme, mass arry, position array, velocity array, acceleration array,
jerk array, snap array, crackle array, timestep array, time array, softening parameter, 
accuracy parameter eta, minimum time (for output), integer specifying whether or not to fix
the timestep
*/
void its_hermite_integrator(int n, int order, double m[], double x[][3], double v[][3], double a[][3],
                        double j[][3], double s[][3], double c[][3], double dt[], double t[], 
                        double eps2, double eta, double *t_min, int fix_dt, int indii[], int *indexes)
{
    /*Looping through particles to find particles with minimum t_i + dt_i*/
    double global_time; /*Variable to hold minimum time*/
    *indexes = 1; /*Counter to keep track used length of array*/
    indii[0] = 0; /*Initializing minimum time as values for particle at index 0*/
    global_time = dt[0] + t[0];

    for(int i = 1; i < n; i++)
    {
        if((dt[i] + t[i]) < global_time)/*New global minimum found*/
        {
            global_time = dt[i] + t[i]; /*Setting new global minimum value*/
            indii[0] = i; /*Restarting array of indexes*/
            *indexes = 1; /*Note that values past indexes will hold onto old/erraneous information*/
        }
        else if((dt[i] + t[i]) == global_time)/*Particle with degenerate time found*/
        {
            indii[*indexes] = i; /*Adding to array of indexes*/
            ++*indexes; /*Increasing used length of array*/

        }
    }


    /*Predicting positions / velocities of all particles at global time*/
    double x_pred[NMAX][DIM], v_pred[NMAX][DIM]; /*Arrays to hold temporary data*/
    double t_factor, t_factor2, t_factor3; /*Temporary variables to save computational time*/
    for(int i = 0; i < n; i++)
    {
        /*Predefining prefactors to reduce computations*/
        t_factor = global_time - t[i]; /* delta t */
        t_factor2 = sqr(t_factor) * 0.5; /* delta t^2 / 2 */
        t_factor3 = t_factor * t_factor2 * k_3_inv; /* delta t^3 / 6 */
        for(int k = 0; k < DIM; k++) /* j is already being used for jerk array */
        {
            x_pred[i][k] = t_factor3*j[i][k] + t_factor2*a[i][k] + t_factor*v[i][k] + x[i][k];
            v_pred[i][k] = t_factor2*j[i][k] + t_factor*a[i][k] + v[i][k];
        }
    }
    
     /*Storing array values for acceleration and jerk in temporary variables*/
    double a0[NMAX][DIM], j0[NMAX][DIM];
    int index;
    for(int i = 0; i < *indexes; i++)
    {
        index = indii[i]; /*Using temporary variable to make my life simple*/
        for(int k = 0; k < DIM; k++)
        {
            a0[index][k] = a[index][k];
            j0[index][k] = j[index][k];
        }
    }
    
    /**************************************Adding number of corrections specified by order***************************/
    double DT, dtI, dtI2, dtI3; /*Will be using at end*/
    double dt3, dt4;
    for(int o = 0; o < order; o++)
    {
        for(int i = 0; i < *indexes; i++)/*Setting acceleration and jerk to 0*/
        {
            index = indii[i]; /*Using temporary variable to make my life simple*/
            for(int k = 0; k < DIM; k++)
            {
                a[index][k] = 0.0;
                j[index][k] = 0.0; /*Need to do this inside of corrections loop*/
            }
        }
        
        /*Finding new acceleration and jerk for degenerate particles*/
        double dr[DIM], dv[DIM];
        double vDotR, R2, rInv, rInv2, rInv3, vDr_3rInv5;
        double mk_rInv3, v5, mk_vDr_3rInv5;

        for(int i = 0; i < *indexes; i++)
        {
            index = indii[i]; /*Using temporary variable to make life easy*/
            for(int k = 0; k < n; k++)
            {
                if(k == index)
                {
                    continue; /*Skipping this iteration as it's looking at self-interaction*/
                }
                else if(m[k] == 0.0)
                {
                    continue; /*skipping any test particles*/
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
                R2 = mag3_sqr(dr) + eps2;
                rInv2 = 1 / R2; /*Important to divide first and then square root to avoid floating point related error*/
                rInv = sqrt(rInv2);
                rInv3 = rInv2 * rInv;
                vDr_3rInv5 = 3.0 * vDotR * rInv3 * rInv2;
                mk_rInv3 = m[k] * rInv3;
                mk_vDr_3rInv5 = m[k] * vDr_3rInv5;

                for(int l = 0; l < 3; l++)
                {
                    a[index][l] +=   mk_rInv3 * dr[l];
                    j[index][l] +=   mk_rInv3 * dv[l] - mk_vDr_3rInv5 * dr[l];

                }
                /*Technically some calculations are being done (but not counted!) twice, specifically between any particles
                noted on the array indii, couldn't think of a clever way to avoid this cleanly right now but could but an
                area for some speed improvement*/
            }
        }

        /*Calculating corrections to predictions of degenerate particles and updating time*/
        for(int i = 0; i < *indexes; i++)
        {
            /*Calculating snap and crackle to use in corrections to predicted location of particle_index*/
            index = indii[i]; /*Using temporary variable to make my life easier*/
            DT = dt[index];/*Not given that all dt's are the same so need to calc prefactors for each particle*/
            dt3 = sqr(DT) * DT;
            dt4 = DT * dt3;
            dtI = 1 / DT;
            dtI2 = sqr(dtI);
            dtI3 = dtI * dtI2;

            t_factor = global_time - t[index]; /* delta t */
            t_factor2 = sqr(t_factor) * 0.5; /* delta t^2 / 2 */
            t_factor3 = t_factor * t_factor2 * k_3_inv; /* delta t^3 / 6 */

            for(int k = 0; k < 3; k++)
            {
                /*Computing higher order derivatives*/
                s[index][k] = (-6 * (a0[index][k] - a[index][k]) - DT*(4 * j0[index][k] + 2*j[index][k]) ) * dtI2; /*Saving values to be used for timestep adjustment*/
                c[index][k] = (12 * (a0[index][k] - a[index][k]) + 6*DT*(j0[index][k] + j[index][k]) ) * dtI3; /*Saving values to be used for timestep adjustment*/

                /*Adjusting predictions*/

                x_pred[index][k] = x[index][k] + DT * (v[index][k] + DT * (0.5 * a0[index][k] + DT * (k_6_inv * j0[index][k] + 
                                    DT * (k_24_inv * s[index][k] + DT * k_120_inv * c[index][k] ) ) ) );
                v_pred[index][k] = v[index][k] + DT * (a0[index][k] + DT * (0.5 * j0[index][k] + DT * (k_6_inv * s[index][k] + 
                                    DT * k_24_inv * c[index][k] ) ) );
            }


        }

    }/***********************************************************************************************************/

    /*****************************Casting predictions into position and velocity arrays*******************************/
    for(int i = 0; i < *indexes; i++)
    {
        index = indii[i];
        for(int k = 0; k < 3; k++)
        {
            x[index][k] = x_pred[index][k];
            v[index][k] = v_pred[index][k];
        }
    }
    /*Updating particle times*/

    for(int i = 0; i < *indexes; i++)
    {
        t[indii[i]] += dt[indii[i]];
    }
    
    /*Updating timesteps of degenerate particles*/
    if(fix_dt == 0) /*Input specifies that we want code to use variable timestep*/
    {
        double accel, jerk, jerk2, snap, snap2, crackle, dtNew;
        for(int i = 0; i < *indexes; i++)
        {
            index = indii[i];

            /*First finding updated value for snap*/
            for(int k = 0; k < 3; k++)
            {
                s[index][k] += c[index][k] * dt[index]; /*no need to adjust crackle value b/c using 3rd order interpolation*/
            }

            /*Finding magnitudes of acceleration, jerk, snap, and crackle*/
            accel = mag3(a[index]);
            jerk2 = mag3_sqr(j[index]);
            jerk = sqrt(jerk2);
            snap2 = mag3_sqr(s[index]);
            snap = sqrt(snap2);
            crackle = mag3(c[index]); 

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
    *t_min = min_arr(n, t);

    return;
}

/*Function to be called prior to first call of its_hermite_integrator, sets initial t, dt, acceleration,  
jerk, snap and crackle values
Inputs: number of particles, mass arry, position array, velocity array, acceleration array,
jerk array, timestep array, time array, softening parameter, initial accuracy parameter eta_s, 
minimum power allowed for timestep, maximum power allowed for timestep, last parameter fixes the dt
step to the inputted parameter if anything aside from 0.0 is inputted
*/
void initialize_parameters(int n, double m[], double x[][3], double v[][3], double a[][3],
                        double j[][3], double s[][3], double c[][3], double dt[], double t[], 
                        double eps2, double eta_s, int min_pow, int max_pow, double dt_fix)
{   
    double dr[3], dv[3];
    double vDotR, R2, rInv, rInv2, rInv3, rInv5, vDr_3rInv5;
    double mi_rInv3, mk_rInv3, mi_vDr_3rInv5, mk_vDr_3rInv5;
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
    for(int i = 0; i < n - 1; i++)
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
            R2 = mag3_sqr(dr) + eps2;
            rInv2 = 1 / R2; /*Important to divide first and then square root to avoid floating point related error*/
            rInv = sqrt(rInv2);
            rInv3 = rInv2 * rInv;
            rInv5 = rInv3 * sqr(rInv);
            vDr_3rInv5 = 3.0 * vDotR * rInv5;
            mi_rInv3 = m[i] * rInv3;
            mk_rInv3 = m[k] * rInv3;
            mi_vDr_3rInv5 = m[i] * vDr_3rInv5;
            mk_vDr_3rInv5 = m[k] * vDr_3rInv5;

            for(int l = 0; l < 3; l++)
            {
                
                a[i][l] +=   mk_rInv3 * dr[l];
                a[k][l] += - mi_rInv3 * dr[l];

                j[i][l] +=   mk_rInv3 * dv[l] - mk_vDr_3rInv5 * dr[l];
                j[k][l] += - mi_rInv3 * dv[l] + mi_vDr_3rInv5 * dr[l];

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
            unbound_dt = eta_s * mag3(a[i]) / mag3(j[i]);
            
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
        set_darr(n, dt, dt_fix);
    }

    return;
}


/*NOTE, WHEN USING FIXED TIMESTEPS, THIS FUNCTION IS UNNEEDED AND SEEMS TO CAUSE HIGHER ERRORS*/
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
    double dt;
    /*Predicting positions / velocities of all particles at global time*/
    for(int i = 0; i < n; i++)
    {
        /*Predefining prefactors to reduce computations*/
        dt =  output_time - t[i];

        for(int k = 0; k < 3; k++) /* j is already being used for jerk array */
        {
            x_pred[i][k] = x[i][k] + dt * (v[i][k] + dt * (0.5 * a[i][k] + dt * (k_6_inv * j[i][k] + 
                            dt * (k_24_inv * s[i][k] + dt * (k_120_inv * c[i][k])))));

            v_pred[i][k] = v[i][k] + dt * (a[i][k] + dt * (0.5 * j[i][k] + dt * (k_6_inv * s[i][k] + 
                            dt * (k_24_inv * c[i][k]))));
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
