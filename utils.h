#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <time.h>
#define NS_PER_SECOND 1e9

#define DATA_OUTPUT 1  
#define DAT_PATH "heat-dat-files"

#define IMAX 400
#define JMAX 400

/* Domain definition */ 
#define TR 50.0   
#define TL 100.0

/* Computation parameters */
#define NMAX 1e6
#define NREP 1       



char* trim(char *s);
void compExactSolution(float **Te, float *x, float kappa, float time);
void DataOutput(  char testname[200], int timestep, float time, float* x, float* y, float **T, float **Te );
void Time_Output( float elapsed_time );
void sub_timespec ( struct timespec t1, struct timespec t2, struct timespec *td);
double simple_sub_timespec ( struct timespec t1 , struct timespec t2 );
void print_tn(float **mat );
void initConditions( float **mat, float *x, float xD);
void allocRows(float **Tn, float** Tn1, float **Te);
