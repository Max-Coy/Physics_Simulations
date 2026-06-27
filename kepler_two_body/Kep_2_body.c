/********************************************
Kep_2_body.c
Max Coy, 7/23
Simple code for testing hermite_integrator.c
*********************************************/

#define NMAX 32
#define DIM 3
#include "hermite_PECn.h"

int main(void)
{    
    /*Initializing variables*/
    double x[NMAX][DIM], v[NMAX][DIM], a[NMAX][DIM], j[NMAX][DIM], s[NMAX][DIM], c[NMAX][DIM]; /*Position and its derivatives*/
    double m[NMAX];
    double output_x[NMAX][DIM], output_v[NMAX][DIM]; /*Temp arrays to be used for outputs*/
    double eta = 0.005, eta_s = 0.001; /*Control accuracy / timestep size */
    double eps2 = pow(2,-10); /*Potential Softening Parameter*/
    double t[NMAX], dt[NMAX]; /*arrays to hold particle times and timesteps respectively*/
    double energy, r_v, t_min;
    double x_mag, v2;
    int indexes;
    int indii[NMAX];
    /*Initializing 2 body system*/
    m[0] = 1.0;
    
    double sm_a, incl, ap, ra, MA, mu; /*Defining orbital parameters to describe initial state*/
    double e[7] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}; /*Testing multiple eccentricities*/ // 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0
    // double e = 0.1;
    sm_a = 1.0; /*Semi-major axis*/
    incl = 0.0; /*inclination*/
    ap = 0.0; /*argument of pericenter*/
    ra = 0.0; /*right ascension of ascending node*/
    MA = k_pi; /*mean anomaly*/
    mu = G *m[0]; /*gravitational parameter*/
    int n = 7; /*number of particles*/
    double tpp;

    for(int i = 0; i < (n-1); i++)/*initializing test particles / leaving first one as center object*/
    {
        m[i+1] = 0.0;
        oe2cart(x[i+1], v[i+1], sm_a, e[i], incl, ap, ra, MA, mu); 
    }

    for(int i = 0; i < 3; i++) /*Initializing central mass*/
    {
        x[0][i] = 0.0;
        v[0][i] = 0.0;
    }

    double period = k_2pi*sm_a * sqrt(sm_a / (m[0] + m[1]));
    
    t_min = 0.0;
    double index = 1.0;
    double resolution = 0.125; ///period/32; /*Outputting data n times per period (not related to timesteps)*/

    double ts = pow(2,-5); /*Fixed timestep*/
    double Time = 1000 * period; /*Total time to be integrated*/
    

    initialize_parameters(n, m, x, v, a, j, s, c, dt, t, eps2, eta_s, 4, 16, ts);
    
    /*Manually varying timesteps*/

    dt[1] = pow(2,-3);
    dt[2] = pow(2,-4);
    dt[3] = pow(2,-5);
    dt[4] = pow(2,-6);
    dt[5] = pow(2,-7);
    dt[6] = pow(2,-8);
    dt[7] = pow(2,-9); /*currently being skipped due to memory constraints*/

    /*Outputting initial state*/

    FILE *stream[7];
    char fn[25];
    for(int i = 0; i < 7; i++)
    {
        sprintf(fn,"timestep_error/ap_dt%d.csv",3+i);
        stream[i] = fopen(fn, "w");
    }
    
    int order = 3;
    int dex;
    /*FULL OUTPUT*/
    double place; /*to be used as a stand in for e in cart2oe() function call*/
    for(int k = 0; k < 7; k++)
    {
        dex = k+1;
        /*Outputing x/v values*/
        fprintf(stream[k], "%1.10f %1.10f %1.0f ", x[dex][0], x[dex][1], x[dex][2]);
        fprintf(stream[k], "%1.0f %1.0f %1.0f ", v[dex][0], v[dex][1], v[dex][2]);
        /*outputing specific energy*/
        x_mag = mag3(x[dex]);
        v2 = sqr(v[dex][0]) + sqr(v[dex][1]) + sqr(v[dex][2]);
        energy = 0.5 * v2 - mu / x_mag;
        fprintf(stream[k], "%1.0f ",energy);
        /*outputing orbital elements*/
        fprintf(stream[k], "%1.10f %1.0f %1.0f %1.12f %1.0f %1.0f", sm_a, e[dex-1], incl, ap, ra, MA);
        fprintf(stream[k], ","); /*Allows me to just comment out outputs i don't want*/
        fprintf(stream[k], "\n");
    }
    print_arr(n, e, 30, 0, 1, 0);
    
    
    double output[NMAX];
    do
    {
    //     // printf("t_min: %1.6f\n",t_min);
        its_hermite_integrator(n, order, m, x, v, a, j, s, c, dt, t, eps2, eta, &t_min, 1,indii, &indexes);

    //     if(0 > 1)//(t_min >= index * resolution)
    //     {
    //         // its_output(output_x, output_v, index*resolution, n, x, v, a, j, s, c, t);
            
    //         /*Full output*/
    //         for(int i = 1; i < n; i++)
    //         {
    //             // // fprintf(stream, "%1.6f %1.6f %1.6f ", output_x[i][0], output_x[i][1], output_x[i][2]);
    //             // // fprintf(stream, "%1.6f %1.6f %1.6f ", output_v[i][0], output_v[i][1], output_v[i][2]);
    //             // // fprintf(stream, "%1.4f %1.4f %1.4f ", x[i][0], x[i][1], x[i][2]);
    //             // // fprintf(stream, "%1.4f %1.4f %1.4f ", v[i][0], v[i][1], v[i][2]);

    //             // /*outputing specific energy*/
    //             // x_mag = mag3(x[i]);
    //             // v2 = sqr(v[i][0]) + sqr(v[i][1]) + sqr(v[i][2]);
    //             // energy = 0.5 * v2 - mu / x_mag;
    //             // // energy = specific_energy(x[i], v[i], mu);
    //             // fprintf(stream, "%1.12f ",energy);
    //             // /*outputing orbital elements*/
    //             // // cart2oe(x[i], v[i], &sm_a, &place, &incl, &ap, &ra, &MA, mu);
    //             // calc_orbital_elements(x[i], v[i], t_min, &sm_a, &place, &incl, &ra, &ap, &tpp);
    //             // e[i-1] = place;
    //             // fprintf(stream, "%1.10f %1.8f %1.2f %1.12f %1.8f %1.10f", sm_a, e[i-1], incl, ap, ra, tpp);
    //             // fprintf(stream,",");
    //         }
    //         // fprintf(stream, "\n");
            
    //         index++;
    //     }
      
        /*Full output*/
        for(int k = 0; k < indexes; k++)
        {
            if(indii[k] == 0)
            {
                continue;
            }
            
            dex = indii[k] ;
            // printf("outputting to %d\t\t",dex);
            fprintf(stream[dex - 1], "%1.10f %1.10f %1.0f ", x[dex][0], x[dex][1], x[dex][2]);
            fprintf(stream[dex - 1], "%1.0f %1.0f %1.0f ", v[dex][0], v[dex][1], v[dex][2]);
            /*outputing specific energy*/
            energy = specific_energy(x[dex], v[dex], mu);
            fprintf(stream[dex - 1], "%1.0f ",energy);
            /*outputing orbital elements*/
            cart2oe(x[dex], v[dex], &sm_a, &place, &incl, &ap, &ra, &tpp, t_min, mu);
            // calc_orbital_elements(x[dex], v[dex], t_min, &sm_a, &place, &incl, &ra, &ap, &tpp);
            e[dex-1] = place;
            fprintf(stream[dex - 1], "%1.10f %1.0f %1.0f %1.12f %1.0f %1.0f", sm_a, e[dex-1], incl, ap-k_2pi, ra, tpp);
            fprintf(stream[dex - 1],",");
            fprintf(stream[dex - 1], "\n");
        }

    }while(t_min < Time);
    

    for(int i = 0; i < 7; i++)
    {
        fclose(stream[i]);
    }
    print_arr(n, e, 30, 0, 1, 1);
    return 0;
}
