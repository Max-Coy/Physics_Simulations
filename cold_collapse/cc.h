/**************************************
cc.h 
Max Coy, 7/23
Header file for cold collapse code
**************************************/

#ifndef NMAX /*Max number of particles allowed*/
#   define NMAX 1024
#endif

#ifndef G /*Gravitational Constant*/
#   define G 1
#endif

#ifndef M /*Total mass of system*/
#   define M 1
#endif

#define drand48() ((double) rand() / RAND_MAX )

/*Importing standard libraries*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

/*Importing Custom Files*/
#include "mfunc.c"          /*Extra math functions*/
#include "leap_frog.c"      /*Leap_frog numerical integrator*/
#include "mkrsphere.c"      /*Generates uniform random spherical distributions*/
