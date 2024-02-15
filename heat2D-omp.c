#include "utils.h"
#include <omp.h>
#include <math.h>
#include <string.h>


int main(int argc, char *argv[]) {

    char test_name[200];
    int i,j,n;

    float xL, xR, yT, yB;                      //domain definition
    float xD;                                  //location of discontinuity
    float dx, dx2, dy, dy2;                    //mesh size and squares
    float x[IMAX], y[IMAX];                    //vertex coords         

    float CFL;                                 //number for stabilty condition(CFL<1)
    float time,dt,tend;                                //current time

    float tio, dtio;                                //output time step

    struct timespec time_start, time_end;              //clock times
    double elapsed, mean, sum;                    //time variables

    float exact;                            //exact solution
    float d,kappa;                             //heat conduction coefficient

    /* 
    COMPUTATION SETTING */

    strcpy(test_name,"Heat2D_explicit");

    xL = 0.0;
    xR = 2.0;
    yB = 0.0;
    yT = 2.0;
    xD = 1.0;

    time = 0.0;
    tend = 0.05;
    tio = 0.0;
    dtio = 0.01;


    CFL = 0.99;
    d = 0.25;
    kappa = 1.0;

    float **Tn, **Tn1, **Te;
    Tn  = allocateMatrix();
    Tn1 = allocateMatrix();
    Te  = allocateMatrix();
    
    
    /* Domain definition on dimension x & y */ 
    dx = (xR - xL) / (float)(IMAX - 1);
    dy = (yT - yB) / (float)(JMAX - 1);

    dx2 = dx * dx;
    dy2 = dy * dy;
    x[0] = xL; y[0] = yB;

    for (i = 0; i < (IMAX - 1); i++)
        x[i + 1] = x[i] + dx;
    for (j = 0; j < (JMAX - 1); j++)
        y[j + 1] = y[j] + dy;

puts("\n  PARALLEL");
puts(" ======================================== ");
#ifdef _USE_OMP
    #pragma omp parallel
    {
        #pragma omp single
        {
        printf(" * Active th : [%d]\t * IMAX: %d\t*\n",omp_get_num_threads(), IMAX);
        printf(" * Cancell   : [%d]\t * JMAX: %d\t*\n",omp_get_cancellation(), JMAX);
        printf(" * \t\t\t * NREP: %d\t*\n", NREP);
        }
    }
#endif
puts(" ======================================== ");
printf(" repetition: "); fflush(stdout);

    int rep;
    for (rep = 0; rep < NREP; rep++) {


        /* Initial condition definition */ 
        initConditions(Tn,Tn1,Te,x,xD, &time);

/*************************************************************************************************/
        clock_gettime(CLOCK_REALTIME, &time_start);
/*************************************************************************************************/
    
/* You can't terminate a parallel construct prematurely. OpenMP has no construct for this and it specifies
 * that parallel regions may have only one exit point (so no branching out of the region...)
 */
#ifdef _USE_OMP
    #pragma omp parallel
    {
#endif
        /* COMPUTATION */
        for (n = 0; n < NMAX; n++) {

            if (time >= tend){ 
                #pragma omp cancel parallel
            }

#ifdef _USE_OMP
    #pragma omp single 
    {
#endif
            /* Computing time step */
            dt = CFL * d * fmax(dx2, dy2) / kappa;
            if ( (time + dt ) > tend) {
                dt = tend - time;
                tio = tend;
            }
        #ifdef _OUTPUT
            if ( ( time + dt ) > tio ) {
                dt = tio - time;
            }
        # endif
            time = time + dt;     // updating time

#ifdef _USE_OMP
    } /* end omp single 1 */
    #pragma omp barrier
#endif

            /* Schema ESPLICITO */
#ifdef _USE_OMP
    #pragma omp for private(i,j,exact) 
#endif
            for (i = 1; i < IMAX - 1; i++) {

                /* computing exact solution for boundaries */ 
                exact = ( (TR + TL) / 2.0)  +  ( erf((x[i] - 1) / (2 * sqrt(kappa * time))) * (TR - TL) / 2.0 );
                
                Tn1[i][0] = exact;
                Tn1[i][JMAX - 1] = exact;
                
                for (j = 1; j < JMAX - 1; j++) {

                    Tn1[i][j] = Tn[i][j] 
                                + ((kappa * dt / dx2) * ( Tn[i + 1][j]  - 2.0 * Tn[i][j] + Tn[i - 1][j]))
                                + ((kappa * dt / dy2) * ( Tn[i][j + 1]  - 2.0 * Tn[i][j] + Tn[i][j - 1]));
                }
                 
            }

#ifdef _USE_OMP
    #pragma omp barrier
    #pragma omp for private(i,j)
#endif
            for (i = 1; i < IMAX - 1; i++) {
                for (j = 0; j < JMAX; j++) {
                    Tn[i][j] = Tn1[i][j];
                }
            }
#ifdef _USE_OMP
    #pragma omp barrier
    #pragma omp single 
    {
#endif

#ifdef _OUTPUT
            if ( (fabs(time - tio) < 1e-12) ) {

                printf("Plotting data output at time %f\n", time);
                /* Computing exact solution in Te */
                compExactSolution(Te, x, kappa, time);
                DataOutput(test_name, n, time, x, y, Tn, Te);

                tio += dtio;
            }
#endif

#ifdef _USE_OMP
    } /*  omp single  */
    #pragma omp barrier
    #pragma omp cancellation point parallel
#endif
        }   /* ==== end for (n) ==== */
#ifdef _USE_OMP
    }
    /* End omp parallel*/
#endif
/*************************************************************************************************/
        clock_gettime(CLOCK_REALTIME, &time_end);
/*************************************************************************************************/

        elapsed = simple_sub_timespec(time_start, time_end);
        sum += elapsed;
        printf("%d  ", rep);
        fflush(stdout);
    
    } /* ==== end for (rep) ==== */

    
    mean = sum / (float)NREP;
    printf("\n\n  >>> Computation time mean :");
    printf( "%.8f\n", mean);
    
    /* Deallocation */ 
    free(Tn); 
    free(Tn1); 
    free(Te);

    if (n == NMAX) {
        printf("Reached maximum number of timesteps\n");
        exit(EXIT_SUCCESS);
    }
    printf("  >>> Exit with no error,experiment succesful \n\n");
    return 0;
}

