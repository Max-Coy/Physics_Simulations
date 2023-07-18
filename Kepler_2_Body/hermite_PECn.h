/*******************************************************
hermite_PECn.h
Max Coy, 7/23
header file for codes using hermite (PEC)^n scheme 
********************************************************/

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

/*Importing standard libraries*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*Importing Custom Files*/
#include "mconst.h"                 /*physical and mathematical consants*/
#include "mfunc.c"                  /*extra math functions*/
#include "array_func.c"             /*Array manipulation and printing functions*/
#include "OE_converter.c"           /*converts between orbital elements and cartesian coordinates*/
#include "hermite_integrator.c"     /*(PEC)^n hermite code*/
