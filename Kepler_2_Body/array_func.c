/****************************************************************************************
array_func.c
Max Coy, 7/23
Some functions for manipulating arrays and printing data
The primary focus for these functions is for setting up initial conditions and debugging
and so clarity will be focused on over speed/efficiency

Functions:
    print_arr()                     print out elements of 1-D array
    print_pt()                      Specialized print function for use with hermite_integrator.c
    print_particles()               Specialized print function for use with hermite_integrator.c

    min_arr()                       returns minimum value in 1-D double array
    min_index()                     returns index of minimum value in 1-D double array
    max_arr()                       returns maximum value in 1-D double array
    max_index()                     returns index of maximum value in 1-D double array

    shift_arr_indicies2()           Shift index location of elements in 2-D array
    shift_arr_values2()             Shift the elements of a 2-D array by fixed amount
    scale_arr_values2()             Scale the elements of a 2-D array by fixed amount

    add_arr()                       add 2 1-D arrays
    add_arr2()                      add 2 2-D arrays
    sub_arr()                       subtract 2 1-D arrays
    sub_arr2()                      subtract 2 2-D arrays
    mult_arr()                      elementwise multiply 2 1-D arrays
    mult_arr2()                     elementwise multiply 2 2-D arrays
    div_arr()                       elementwise divide 2 1-D arrays
    div_arr2()                      elementwise divide 2 2-D arrays

    cast_iarr_darr()                Convert 1-D int array to double array
    cast_darr_iarr()                Convert 1-D double array to int array
    set_darr()                      Set all values in double array
    set_iarr()                      Set all values in int array

*****************************************************************************************/

/*Prints out 1D array horizontally
Inputs: length of array, array, integer specifying precision, integer specifying print type, 
integer specifying direction, int specifying horizontal buffer (only for vertical printing)

    prec is allowed to range from 0 to 99

    type = 0 > printing in float
    type = 1 > printing as SN

    dir = 0 > printing horizontally
    dir = 1 > printing vertically

    buff is allowed to range from 0 to 9, it adds \t buff times prior to printing
*/
void print_arr(int len, double arr[], int prec, int type, int dir, int buff)
{
    char out[7] = "%1."; /*Will be concatenated with places*/
    char places[4]; /*First two slots are for precision, last is for type*/

    /*Making sure specified precision is within [0,99] to prevent index out of bounds*/
    if(prec > 99) /*Precision is too large*/
    {
        printf("\n\nWarning, specified precision is to high, assigning 99 instead\n\n");
        sprintf(places,"%d",99);
    }else if(prec < 0)/*Precision is too small*/
    {
        printf("\n\nWarning, specified precision is too small, assigning 0 instead\n\n");
        sprintf(places,"%d",0);
    }else{/*Precision is within bounds, adding to output string*/
        sprintf(places,"%d",prec);
    }

    /*Setting up output string type*/
    if(type == 0)/*Printing as float*/
    {
        strcat(places,"f");
    }else if(type == 1)/*Printing in SN*/
    {
        strcat(places,"e");
    }else{/*Invalid type specified, setting as SN*/
        printf("\n\nWarning specified print type is not valid, setting as SN\n\n");
        strcat(places,"e");
    }
    
    strcat(out,places); /*Forming output string*/

    /*Printing out data*/
    if(dir == 0) /*Printing Data Horizontally*/
    {
        for(int i = 0; i < len; i++)
        {
            printf(out,arr[i]); /*printing data entry*/
            printf("\t"); /*Adding horizontal space between entries*/
        }
        return;
    }
    
    if(buff < 0)/*Checking buffer values*/
    {
        printf("\n\nWarning, specified buffer is too small, setting to 0\n\n");
        buff = 0;
    }else if(buff > 9)
    {
        printf("\n\nWarning, specified buffer is too large, setting to 9\n\n");
        buff = 9;
    }
    
    if(dir == 1) /*Printing Data Vertically*/
    {
        for(int i = 0; i < len; i++)
        {
            for(int j = 0; j < buff; j++)/*Adding horizontal buffer*/
            {
                printf("\t"); 
            }
            printf(out,arr[i]);  /*printing data entry*/
            printf("\n"); /*Adding vertical space between entries*/
        }
    }else{/*Invalid direction given, printing vertically*/
        printf("\n\nWarning specified print direction is not valid, printing vertically\n\n");
        for(int i = 0; i < len; i++)
        {
            for(int j = 0; j < buff; j++)/*Adding horizontal buffer*/
            {
                printf("\t"); 
            }
            printf(out,arr[i]);  /*printing data entry*/
            printf("\n"); /*Adding vertical space between entries*/
        }
    }
    
    return;
}

/*Prints out particle data
Inputs: particle index, arrays position and all of its derivatives in ascending order, number of dimensions, 
desired display precision, display type, horizontal display buffer (for vertical printing only)

    prec is allowed to range from 0 to 99

    type = 0 > printing in float
    type = 1 > printing as SN

    dir = 0 > printing horizontally
    dir = 1 > printing vertically

    buff is allowed to range from 0 to 9, it adds \t buff times prior to printing
*/
void print_pt(int index, double x[], double v[], double a[], double j[], double s[], double c[], 
              int dim, int prec, int type, int buffer)
{
    int direction = 1; /*Printing vertically*/
    printf("Particle: %1.0d\n",index);
    printf("\tX:\n");
    print_arr(dim,x,prec,type,direction,buffer);
    printf("\tV:\n");
    print_arr(dim,v,prec,type,direction,buffer);
    printf("\tA:\n");
    print_arr(dim,a,prec,type,direction,buffer);
    printf("\tJ:\n");
    print_arr(dim,j,prec,type,direction,buffer);
    printf("\tS:\n");
    print_arr(dim,s,prec,type,direction,buffer);
    printf("\tC:\n");
    print_arr(dim,c,prec,type,direction,buffer);
    printf("\n");
    return;
}

/*Prints out 3 dimensional particle data for list of particles
Inputs: number of particles, arrays of particle position and all of its derivatives in ascending order,
number of dimensions desired display precision, display type, horizontal display buffer (for vertical printing only)

    prec is allowed to range from 0 to 99

    type = 0 > printing in float
    type = 1 > printing as SN

    dir = 0 > printing horizontally
    dir = 1 > printing vertically

    buff is allowed to range from 0 to 9, it adds \t buff times prior to printing
*/
void print_particles(int n, double x[][DIM], double v[][DIM], double a[][DIM], double j[][DIM], 
                     double s[][DIM], double c[][DIM], int dim, int precision, int type, int buffer)
{
    
    for(int i = 0; i < n; i++)
    {
        print_pt(i, x[i], v[i], a[i], j[i], s[i], c[i], dim, precision, type, buffer);
        printf("\n");
    }
    return;
}


/*Returns minimum value of 1D double array input given n its length*/
double min_arr(int n, double input[])
{
    double min = input[0];
    for(int i = 1; i < n; i++)
    {
        min = fmin(min, input[i]);
    }
    return min;
}

/*Returns index of minimum value of 1D double array input given n its length*/
int min_index(int n, double input[])
{
    double min;
    int index;
    min = input[0];
    index = 0;
    for(int i = 1; i < n; i++)
    {
        if( fmin(min, input[i]) != min)
        {
            index = i;
            min = input[i];
        }
    }
    return index;
}

/*Returns maximum value of 1D double array input given n its length*/
double max_arr(int n, double input[])
{
    double max = input[0];
    for(int i = 1; i < n; i++)
    {
        max = fmax(max, input[i]);
    }
    return max;
}

/*Returns index of minimum value of 1D double array input given n its length*/
int max_index(int n, double input[])
{
    double max;
    int index;
    max = input[0];
    index = 0;
    for(int i = 1; i < n; i++)
    {
        if( fmax(max, input[i]) != max)
        {
            index = i;
            max = input[i];
        }
    }
    return index;
}

/*Adjust position of informations in 1D array
Inputs: number of particles, Current offset, New offset, array to adjust, number of dimensions
Warning: poor indexing can lead to data being overwritten
*/
void shift_arr_indicies(int n, int current_offset, int new_offset, double arr[])
{
    /*Ensuring there is space in the array for new particles*/
    if(n + new_offset > NMAX)
    {
        printf("\n\nWarning, attempting to index variables out of array bounds, terminating call\n\n");
        return;
    }
    /*Adjusting information*/
    for(int i = 0; i < n; i++)
    {           
        arr[i+new_offset] = arr[i+current_offset];
        arr[i+current_offset] = 0;
    }
    return;
}

/*Adjust position of informations in 2D array
Inputs: number of particles, Current offset, New offset, array to adjust, number of dimensions
Warning: poor indexing can lead to data being overwritten
*/
void shift_arr_indicies2(int n, int current_offset, int new_offset, double arr[][DIM], int dim)
{
    /*Ensuring there is space in the array for new particles*/
    if(n + new_offset > NMAX)
    {
        printf("\n\nWarning, attempting to index variables out of array bounds, terminating call\n\n");
        return;
    }
    /*Adjusting information*/
    for(int i = 0; i < n; i++)
    {           
        for(int j = 0; j < dim; j++)
        {
            arr[i+new_offset][j] = arr[i+current_offset][j];
            arr[i+current_offset][j] = 0;
        } 
    }
    return;
}

/*Shift the specified values in 2D array by fixed amount
Inputs: number of particles, index offset, array containing desired offset,
array to be modified, number of dimensions
*/
void shift_arr_values2(int n, int offset, double shift[DIM], double arr[][DIM], int dim)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            arr[i + offset][j] += shift[j];
        }
    }
    return;
}

/*Scale the specified values in 2D array by fixed amount 
Inputs: number of particles, index offset, array containing desired scaling,
array to be modified, number of dimensions
*/
void scale_arr_values2(int n, int offset, double scaling[DIM], double arr[][DIM], int dim)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            arr[i + offset][j] *= scaling[j];
        }
    }
    return;
}

/*Add two 1D arrays together
Takes length of arrays, arrays to be added, output array
*/
void add_arr(int n, double arr1[], double arr2[], double out[])
{
    for(int i = 0; i < n; i++)
    {
        out[i] = arr1[i] + arr2[i];
    }
    return;
}

/*Add two 2D arrays together
takes length of arrays, arrays to be added, output array, dimension of 2nd axis
*/
void add_arr2(int n, double arr1[][DIM], double arr2[][DIM], double out[][DIM], int dim)
{
    for(int i = 0; i < n; i++)
    {
        add_arr(dim, arr1[i], arr2[i], out[i]);
    }
    return;
}

/*Subtract two 1D arrays
Takes length of arrays, arrays to be subtracted, output array
*/
void sub_arr(int n, double arr1[], double arr2[], double out[])
{
    for(int i = 0; i < n; i++)
    {
        out[i] = arr1[i] - arr2[i];
    }
    return;
}

/*Subtract two 2D arrays 
takes length of arrays, arrays to be subtracted, output array, dimension of 2nd axis
*/
void sub_arr2(int n, double arr1[][DIM], double arr2[][DIM], double out[][DIM], int dim)
{
    for(int i = 0; i < n; i++)
    {
        sub_arr(dim, arr1[i], arr2[i], out[i]);
    }
    return;
}

/*Multiply two 1D arrays together
Takes length of arrays, arrays to be multiplied, output array
*/
void mult_arr(int n, double arr1[], double arr2[], double out[])
{
    for(int i = 0; i < n; i++)
    {
        out[i] = arr1[i] * arr2[i];
    }
    return;
}

/*Multiply two 2D arrays together
takes length of arrays, arrays to be multiplied, output array, dimension of 2nd axis
*/
void mult_arr2(int n, double arr1[][DIM], double arr2[][DIM], double out[][DIM], int dim)
{
    for(int i = 0; i < n; i++)
    {
        mult_arr(dim, arr1[i], arr2[i], out[i]);
    }
    return;
}

/*Divide two 1D arrays together
Takes length of arrays, arrays to be divided, output array
*/
void div_arr(int n, double arr1[], double arr2[], double out[])
{
    for(int i = 0; i < n; i++)
    {
        out[i] = arr1[i] / arr2[i];
    }
    return;
}

/*Divide two 2D arrays together
takes length of arrays, arrays to be divided, output array, dimension of 2nd axis
*/
void div_arr2(int n, double arr1[][DIM], double arr2[][DIM], double out[][DIM], int dim)
{
    for(int i = 0; i < n; i++)
    {
        div_arr(dim, arr1[i], arr2[i], out[i]);
    }
    return;
}

/*Casts an array of int values into an array of double values*/
void cast_iarr_darr(int n, int arr1[], double arr2[])
{
    for(int i = 0; i < n; i++)
    {
        arr2[i] = (double) arr1[i];
    }
    return;
}

/*Casts an array of double values into an array of int values*/
void cast_darr_iarr(int n, double arr1[], int arr2[])
{
    for(int i = 0; i < n; i++)
    {
        arr2[i] = (int) arr1[i];
    }
    return;
}

/*Set every value in double array*/
void set_darr(int n, double arr[], double val)
{
    for(int i = 0; i < n; i++)
    {
        arr[i] = val;
    }
    return;
}

/*Set every value in int array*/
void set_iarr(int n, int arr[], int val)
{
    for(int i = 0; i < n; i++)
    {
        arr[i] = val;
    }
    return;
}
