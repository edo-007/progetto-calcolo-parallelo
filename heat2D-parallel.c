#include "utils.h"
#include <omp.h>
#include <math.h>
/* 
    ^
    | __________________  
    |[_______50_________] 
    |[ |              | ]
    |[e|              |e]
    |[ |              | ]
    |[_|______________| ]
  __|[_______100________]____>
    |
 */
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
    int plotindex;

    struct timespec time_start, time_end;              //clock times
    double elapsed, mean, sum;                    //time variables

    float Te_point;                            //exact solution
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

    plotindex = 0;

    CFL = 0.99;
    d = 0.25;
    kappa = 1.0;
    
    puts("\n  ** Parallel *****************");
    
    float **Tn, **Tn1, **Te;
    Tn  = (float**)malloc(IMAX * sizeof(float*));
    Tn1 = (float**)malloc(IMAX * sizeof(float*));
    Te  = (float**)malloc(IMAX * sizeof(float*));
    allocRows(Tn, Tn1, Te);
    
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

    int rep;
    for (rep = 0; rep < NREP; rep++) {
    
        /* Initial condition definition */ 
        initConditions(Tn,x,xD);
        clock_gettime(CLOCK_REALTIME, &time_start);


#ifdef _USE_OMP
    #pragma omp parallel private(i, j, Te) shared(time)
    {
        // #pragma omp single 
        // {
        //     printf("   - number th : %d\n", omp_get_num_threads());
        //     printf("   - max th : %d\n\n", omp_get_max_threads());
        // }
#endif
        /* COMPUTATION */
        for (n = 0; n < NMAX; n++) {                // nnnnnnnnnnnnnnnnnnnnnnnnnnnnn
            if (time >= tend) {
                printf("%f\n", time);
                break;
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
    #pragma omp for 
#endif
            for (i = 1; i < IMAX - 1; i++) {

                Te_point = ( (TR + TL) / 2.0)  +  ( erf((x[i] - 1) / (2 * sqrt(kappa * time))) * (TR - TL) / 2.0 );
                
                Tn1[i][0] = Te_point;
                Tn1[i][JMAX - 1] = Te_point;
                
                for (j = 1; j < JMAX - 1; j++) {

                    Tn1[i][j] = Tn[i][j] 
                                + ((kappa * dt / dx2) * ( Tn[i + 1][j]  - 2.0 * Tn[i][j] + Tn[i - 1][j]))
                                + ((kappa * dt / dy2) * ( Tn[i][j + 1]  - 2.0 * Tn[i][j] + Tn[i][j - 1]));
                }
                 /* computing exact solution for boundaries */ 
            }

#ifdef _USE_OMP
    #pragma omp barrier
    /* Overwrite */
    #pragma omp for
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
#endif
        }   /* ==== end for (n) ==== */
#ifdef _USE_OMP
        printf("thread %d\n", omp_get_thread_num());
    }
    /* End omp parallel*/
#endif
        clock_gettime(CLOCK_REALTIME, &time_end);
        elapsed = simple_sub_timespec(time_start, time_end);
        sum += elapsed;
    
    } /* ==== end for (rep) ==== */

    puts("  ====================");


    /* Output */
    mean = sum / (float)NREP;

    printf("Final time output: %.8f\n", mean);
    // Time_Output(mean);

    /* Deallocation */ 
    free(Tn);
    free(Tn1);
    free(Te);

    if (n == NMAX) {
        printf("Reached maximum number of timesteps\n");
        exit(EXIT_SUCCESS);
    }
    printf("Exit with no error,experiment succesful \n");
    return 0;
}

