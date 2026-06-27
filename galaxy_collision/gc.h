/*************************************
gc.h
Max Coy, 7/23
header file for galaxy collison code
*************************************/

#ifndef NMAX /*Max number of particles allowed*/
#   define NMAX 16 
#endif

#ifndef G /*Gravitational Constant*/
#   define G 1
#endif

#ifndef M /*Total mass of system*/
#   define M 1
#endif

#ifndef DIM /*Number of dimensions*/
#   define DIM 3
#endif

#define drand48() ((double) rand() / RAND_MAX )

/*Importing standard libraries*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

/*Importing Custom Files*/
#include "mconst.h"         /*physical and mathermatical constants*/
#include "mfunc.c"          /*Extra math functions*/
#include "array_func.c"     /*Array manipulation and printing functions*/
#include "leap_frog.c"      /*Leap_frog numerical integrator*/
#include "mkplummer.c"      /*Generates uniform random spherical distributions*/