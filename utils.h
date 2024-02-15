#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#define NS_PER_SECOND 1e9

#define DATA_OUTPUT 1  
#define DAT_PATH "heat-dat-files"

#define IMAX 8
#define JMAX 8

/* Domain definition */ 
#define TR 50.0   
#define TL 100.0

/* Computation parameters */
#define NMAX 1e6
#define NREP 10  

#define BIDIMENSIONAL 2
#define NO_REORDER 0



char* trim(char *s);
void compExactSolution(float **Te, float *x, float kappa, float time);
void DataOutput(  char testname[200], int timestep, float time, float* x, float* y, float **T, float **Te );
void Time_Output( float elapsed_time );
void sub_timespec ( struct timespec t1, struct timespec t2, struct timespec *td);
double simple_sub_timespec ( struct timespec t1 , struct timespec t2 );
void print_tn(float **mat );
void initConditions( float **Tn,float **Tn1, float **Te, float *x, float xD,float *time);
void MPI_Domain_Check(int my_size);
float **allocateMatrix();
float **allocateMatrix_MPI(int total_x, int total_y);



#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */