#include "utils.h"



/* 
    ^
    | __________________  
    |[ |              | ] 
    |[ |     50       | ]
    |[e|______________|e]
    |[ |              | ]
    |[_|     100      | ]
  __|[_|_____ ________|_]__>
    |
 */
int main(int argc, char *argv[]) {

    char test_name[200];
    int i,j,n;

    float xL, xR, yT, yB;                      //domain definition
    float xD;                                  //location of discontinuity
    float dx, dx2, dy, dy2;                    //mesh size and squares
    float x[IMAX], y[JMAX];                    //vertex coords         

    float CFL;                                 //number for stabilty condition(CFL<1)
    float time,dt,tend;                                //current time

    float tio, dtio;                                //output time step

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

    CFL = 0.99;
    d = 0.25;
    kappa = 1.0;


    puts(" Serial ");

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

    int rep;
    for (rep = 0; rep < NREP; rep++) {
    
        /* Initial condition definition */ 
        initConditions(Tn,Tn1,Te,x,xD,&time);
/*************************************************************************************************/
        clock_gettime(CLOCK_REALTIME, &time_start);
/*************************************************************************************************/

        /* COMPUTATION */
        for (n = 0; n < NMAX; n++) {                // nnnnnnnnnnnnnnnnnnnnnnnnnnnnn

            if (time >= tend)
                break;

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

            /* Schema ESPLICITO */
            for (i = 1; i < IMAX - 1; i++) {
                for (j = 1; j < JMAX - 1; j++) {

                    Tn1[i][j] = Tn[i][j] 
                                + ((kappa * dt / dx2) * ( Tn[i + 1][j]  - 2.0 * Tn[i][j] + Tn[i - 1][j]))
                                + ((kappa * dt / dy2) * ( Tn[i][j + 1]  - 2.0 * Tn[i][j] + Tn[i][j - 1]));
                }

                 /* computing exact solution for boundaries */ 
                Te_point = ( (TR + TL) / 2.0)  +  ( erf((x[i] - 1) / (2 * sqrt(kappa * time))) * (TR - TL) / 2.0 );

                Tn1[i][0] = Te_point;
                Tn1[i][JMAX - 1] = Te_point;
            }      
            time = time + dt;     // updating time

            /* Overwrite last solution */
            for (i = 1; i < IMAX - 1; i++) {
                for (j = 0; j < JMAX; j++) {
                    Tn[i][j] = Tn1[i][j];
                }
            }
#ifdef _OUTPUT
            if ( (fabs(time - tio) < 1e-12) ) {

                printf("Plotting data output at time %f\n", time);
                /* Computing exact solution in Te */
                compExactSolution(Te, x, kappa, time);
                DataOutput(test_name, n, time, x, y, Tn, Te);

                tio += dtio;
            }
#endif

        } /* n */ 

/*************************************************************************************************/
        clock_gettime(CLOCK_REALTIME, &time_end);
/*************************************************************************************************/

        elapsed = simple_sub_timespec(time_start, time_end);
        sum += elapsed;
    
    } /* ==== end for (rep) ==== */

    /* Output */
    mean = sum / (float)NREP;
    printf("Final time output: %f\n", mean);

    /* Deallocation */ 
    free(Tn);
    free(Tn1);
    free(Te);
    printf("n = %d\n",n);
    if (n == NMAX) {
        printf("Reached maximum number of timesteps\n");
        exit(EXIT_SUCCESS);
    }
    printf("Exit with no error,experiment succesful \n");
    return 0;
}

