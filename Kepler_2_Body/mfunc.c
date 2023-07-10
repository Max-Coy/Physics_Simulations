/**********************************************************
 mfunc.c
Max Coy, 7/23
Some useful math functions not found in core library
**********************************************************/

#define sqr(x) ((x) * (x)) /*Simple squaring function*/
#define mag3_sqr(x) (sqr(x[0]) + sqr(x[1]) + sqr(x[2])) /*Computes magnitude^2 of 3d vector*/
#define mag3(x) (sqrt(mag3_sqr(x))) /*computes magnitude of 3d vector*/
#define dot3(x,y) (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]) /*computes dot product of 2 3d vectors*/
#define dist3(x,y) (sqrt( sqr(x[0] - y[0]) + sqr(x[1] - y[1]) + sqr(x[2] - y[2]) )) /*computes distance between 2 3d vectors*/

/* Takes 2 3D vectors A and B and saves result to output vector*/
void cross_product(double A[3], double B[3], double output[3])
{
    output[0] = A[1]*B[2] - A[2]*B[1];
    output[1] = A[2]*B[0] - A[0]*B[2];
    output[2] = A[0]*B[1] - A[1]*B[0];
    return;
}
