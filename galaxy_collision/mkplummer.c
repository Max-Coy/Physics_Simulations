/*************************************************************************************************************************
mkplummer.c
Max Coy, 7/23
Generates a plummer sphere distribution in cartesian coordinates

Functions
    calc_mass_cutoff()                      Int         Calculates the mass enclosed by system given cutoff radius
    radial_distrubtion()                    Int         Creates radial distribution of particles
    generate_particle_distribution()        Int         Creates initial particle distribution in cartesian coordinates
    generate_velocity_distribution()        Int         Creates initial velocity distribution
    mkplummer()                             Ext         Generates a Plummer Sphere

External Variables:   
    double  m[]                             O           mass of particles
    int     M                               SC          Total mass of system
    int     n                               I           number of particles in system
    int     NMAX                            SC          max number of particles allowed in simulation
    double  r_cut                           I           cutoff radius for system
    double  v[][]                           O           velocity of particles
    double  x[][]                           O           position of particles

Requires symbolic constants from "mconst.h" to run
WARNING: Currently generate_velocity_distrubition has a fudge factor of 7 in velocity magnitude calculation to prevent
the immeadiate collapse of the system. Not sure where the bug is in the calculation that neccessitates this. 
*************************************************************************************************************************/



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
    double scale_factor = 5.333 * k_pi_inv; /* 16 / 3 pi*/ 
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
    double exp = -2.0/3.0; /*Exponent*/
    for(int i = 0; i < n; i++)
    {
        r[i] = sqrt(pow(drand48()*mass_cut,exp) - 1); 
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
        theta = drand48() * 2*k_pi;
        phi = acos(1 - drand48()*2);
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
    double v1[3], v2[3], vp[2];

    for(int i = 0; i < n; i++)
    {
        do
        {
            y = drand48() * 0.1;
            weight = drand48();
            weight = pow(1-sqr(weight),3.5) * sqr(weight);
        }while(y > weight); /*Half-understand / half don't understand what's happening in the loop lmao*/

        v_esc = k_sqrt2 * pow(1 + sqr(r[i]), -0.25); /*Escape velocity at given radius*/
        v_mag = v_esc*weight * 7; /*Total magnitude of particle velocity*/
        /*^ *7 is a fudge factor added to prevent collapse, not really sure what's going wrong with
        my velocity calc that is making me need it but alas...*/

        /*Pointing velocity in random direction*/
        theta = drand48() * 2*k_pi;
        phi = acos(1 - drand48()*2);

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
