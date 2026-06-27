/*********************************************************************************************
OE_converter.c
Max Coy, 7/23
Transforms orbital elements into cartesian coordinants and vice versa

Functions: 
    specific_energy()                   Calculates the specific energy of particle
    atan2p()                            standard atan2() function but outputs values in [0,2pi]
    MA2EA()                             Converts mean anomaly to eccentric anomaly
    oe2cart()                           Converts orbital elements to cartesian coordinates
    cart2oe()                           Converts cartesian coordinates to orbital elements

External Variables:
    a:                      I/O         semi-major axis
    ap:                     I/O         argument of pericenter
    e:                      I/O         eccentricity
    E:                      O           specific energy
    incl:                   I/O         inclination of orbit
    MA:                     I           mean anomaly
    mu:                     I           gravitational parameter
    ra:                     I/O         right ascension of ascending node
    t:                      I           current time
    tp:                     O           time at pericenter

Requires functions from "mfunc.c" and "array_func.c" to run
Requires constants from "mconst.h" to run
************************************************************************************************/

/*Calculates specific energy, assumes 3d case*/
#define specific_energy(x,v,mu) (mag3_sqr(v) * 0.5 - mu / mag3(x)) 

/*atan2() function but will only return positive radian values*/
double atan2p(double x, double y)
{
    double theta;
    theta = atan2(x,y);
    theta += k_2pi * (theta < 0.0);
    return theta;
}

/*Converts mean anomaly to eccentric anomaly by solving the equation:
            MA = EA - e sin( EA ) 
numerically using halley's method,
Inputs: mean anomaly, eccentricity
*/
double MA2EA(double MA, double e)
{
    double f, df, EN_0; /*temporary variables to reduce total required calculations*/
    double esinE, ecosE;
    double EN_1 = MA; /*Seeding initial guess as the Mean anomaly*/
    double del = 1e-10; /*Window of Convergence*/

    do
    {
        EN_0 = EN_1; /*Saving current guess*/
        esinE = e*sin(EN_0);
        ecosE = e*cos(EN_0);
        f = EN_0 - esinE - MA;
        df = 1 - ecosE;
        EN_1 = EN_0 - (2*f*df) / (2*sqr(df) - f*esinE); /*Calculating new guess*/
        
    }while ( fabs(EN_0-EN_1) > del); /*Checking to see if difference is within allowed window*/

    return EN_1;
}

/*Converts orbital elements to cartesian coordinates
Inputs: arrays to hold position and velocity outputs, semi-major axis, eccentricity, inclination, 
argument of pericenter, right ascension of ascending node, mean anomaly, gravitational parameter
*/
void oe2cart(double x[3], double v[3], double a, double e, double incl, double ap, double ra, 
             double MA, double mu)
{
    double EA = MA2EA(MA,e); /*Eccentric anomaly*/

    double TA = 2 * atan( sqrt( (1+e) / (1-e) ) * tan(EA * 0.5) ); /*True Ascension*/

    double r = a * (1 - e * cos(EA) ); /*radius*/

    double p = a * (1 - sqr(e)); /*Will use this value multiple times*/

    double h = sqrt(mu * p); /*angular momentum*/

    /*Predefining variables to reduce calculations*/
    double TAap, sinTA, sinTAap, cosTAap, sinra, cosra, cosincl, sinincl;
    TAap = TA+ap;
    sinTA = sin(TA);
    sinTAap = sin(TAap);
    cosTAap = cos(TAap);
    sinra = sin(ra);
    cosra = cos(ra);
    cosincl = cos(incl);
    sinincl = sin(incl);

    /*Calculating position vector*/
    x[0] = r * (cosra * cosTAap - sinra * sinTAap * cosincl);
    x[1] = r * (sinra * cosTAap + cosra * sinTAap * cosincl);
    x[2] = r * (sinincl * sinTAap);

    /*Calculating velocity vector*/
    double he_rp, h_r;
    he_rp = h*e / (r*p);
    h_r = h / r;
    v[0] = x[0] * he_rp * sinTA - h_r * (cosra * sinTAap + sinra * cosTAap * cosincl);
    v[1] = x[1] * he_rp * sinTA - h_r * (sinra * sinTAap - cosra * cosTAap * cosincl);
    v[2] = x[2] * he_rp * sinTA + h_r * sinincl * cosTAap;

    return;
}

/*Converts caresian coordinates to orbital elements
Inputs: arrays containing position and velocity values, pointers to variables for semi-major axis, 
eccentricity, inclination, argument of pericenter, right ascension of ascending node, tp, t, and
a variable for the gravitational parameter
Warning: Fails for eccentricities >= 1 (i.e. must be a bound orbit)
*/
void cart2oe(double x[3], double v[3], double *a, double *e, double *incl, double *ap, double *ra,
             double *tp, double t, double mu)
{
    double h[3]; /*Angular momentum*/ 
    cross_product(x,v,h);

    double h_mag = mag3(h); /*Magnitude of angular momentum*/

    double r = mag3(x); /*radius*/
    
    double v2 = sqr(v[0]) + sqr(v[1]) + sqr(v[2]);  /*velocity magnitude^2*/

    double E = v2 * 0.5 - mu / r; /*Specific energy*/


    double hrh[3]; 
    cross_product(v,h,hrh);
    double r_inv = 1 / r;
    double unit_x[3] = {r_inv, r_inv, r_inv};
    mult_arr(3, x, unit_x, unit_x);
    sub_arr(3, hrh, unit_x, hrh);


    *a = - mu / (2*E); /*semi-major axis*/

    *e = sqrt( 1 - sqr(h_mag) / (*a * mu) ); /*eccentricity*/

    *incl = acos(h[2] / h_mag); /*inclination */

    

    double sinincl = sin(*incl);

    if(*incl) /*checking to see if inclination is not 0*/
    {
        *ra = atan2p( h[0], -h[1] ); /*right ascension of ascending node*/   
    }
    else /*Inclination is 0*/
    {
        *ra = 0.0;
    }
    if(*e) /*Checking to see eccentricity is not 0*/
    {
        if(*incl) /*ensuring inclination is nonzero to prevent divide by 0*/
        {
            *ap = atan2p(hrh[2] / sin(*incl), hrh[0] * cos(*ra) + hrh[1] * sin(*ra));
        }
        else /*inclination is 0, using alternate formula*/
        {
            *ap = atan2p(hrh[1], hrh[0]);
        }

        double sqrta = sqrt(*a);
        double esinea = dot3(x,v) / sqrta;
        double EA = atan2p(esinea, 1.0 - r/(*a)); /*eccentric anomaly*/
        *tp = fmod(t - (EA - esinea) * (*a) *sqrta, k_2pi);
    }
    else /*Eccentricity is 0*/
    {
        *ap = 0.0;
        *tp = 0.0;
    }
    
    return;
}
