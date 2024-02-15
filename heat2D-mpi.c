#include "utils.h"
#include <mpi.h>
// #define NDIM 2 
// #define TRUE 1
// #define MPI_ERROR 2

void MPI_Domain_Check(int my_size){

    if ( (my_size % 2) != 0 ){

            printf("Error: number of cpu must be EVEN!\n");

            MPI_Finalize();
            exit(0);
        }
    
        /* (2) Create a cartesian topology */

    if ( (IMAX % 2) != 0 ){
        printf("Error: number of x-cells must be EVEN!");
        MPI_Finalize();
        exit(0);
    }
    if ( (JMAX % 2) != 0 ){
        printf("Error: number of y-cells must be EVEN!");
        MPI_Finalize();
        exit(0);
    }
}


void processToMap(int *xs, int *ys, int *xe, int *ye, int xcell, int ycell, int x_domains, int y_domains) {

   /* Index variables */
   int i, j;

   for (i=0;i<x_domains;i++) {
      ys[i] = 2;
      /* Here, ye(0:(x_domains-1)) = 2+ycell-1 */
      ye[i] = ys[i]+ycell-1;
   }

   for (i=1;i<y_domains;i++)
      for (j=0;j<x_domains;j++) {
         ys[i*x_domains+j] = ys[(i-1)*x_domains+j]+ycell+2;
         ye[i*x_domains+j] = ys[i*x_domains+j]+ycell-1;
      }

   for (i=0;i<y_domains;i++) {
      xs[i*x_domains] = 2;
      xe[i*x_domains] = xs[i*x_domains]+xcell-1;
   }

   for (i=1;i<=y_domains;i++)
      for (j=1;j<x_domains;j++) {
         xs[(i-1)*x_domains+j] = xs[(i-1)*x_domains+(j-1)]+xcell+2;
         xe[(i-1)*x_domains+j] = xs[(i-1)*x_domains+j]+xcell-1;
      }
}
void initConditionsMPI(float **Tn, float **Tn1, int total_x_size, int total_y_size, float *x, float xD, float *time){


   int i, j;

    for (i = 0; i < total_x_size; i++){
        for (j = 0; j < total_y_size; j++){
                if ( x[i] < xD )
                    Tn[i][j] = TL;
                else
                    Tn[i][j] = TR;    
                
                Tn1[i][j] = 0.0;
        }
    }
 }




int main(int argc, char *argv[]) {

    char test_name[200];
    int i,j,n;

    float xL, xR, yT, yB;                      // domain definition
    float xD;                                  // location of discontinuity
    float dx, dx2, dy, dy2;                    // mesh size and squares
    float x[IMAX], y[JMAX];                    // vertex coords         

    float CFL;                                 // number for stabilty condition(CFL<1)
    float time,dt,tend;                        // current time

    float tio, dtio;                           tio = 0; dtio = dtio+tio+1;// output time step 

    struct timespec time_start, time_end;      // clock times
    double elapsed, mean, sum;                 // time variables

    float Te_point;                            // exact solution
    float d,kappa;                             // heat conduction coefficient

    int total_x_size, total_y_size;

    /* MPI  variable  ------------------------------------------------- */
    
    int   my_rank;
    int   my_size;

    int   x_thread;
    int   y_thread;

    int    ndim;
    int    dim   [2] = {0}; 
    int    period[2] = {0}; 
    int    coord [2];
    
    MPI_Status  status;
    MPI_Comm    new_comm;
    MPI_Datatype column_type;
    
    int TCPU, BCPU, RCPU, LCPU;  // neighbor ranks of my_rank


    

    /* end MPI variable ---------------------------------------------------- */


    /* COMPUTATION SETTING */
    strcpy(test_name,"Heat2D_explicit");

    xL = 0.0; xR = 2.0;
    yB = 0.0; yT = 2.0;
    xD = 1.0;

    time = 0.0;
    tend = 0.05;
    tio  = 0.0;
    dtio = 0.01;

    CFL = 0.99; d = 0.25; kappa = 1.0;

    /* (1) MPI inizialization */ 

    if ( MPI_Init(&argc, &argv) != MPI_SUCCESS ) {
        fprintf(stderr , "Error in MPI_Init.\n" );
        exit(-1);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &my_size);       // Size
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
    /* Check for even number of CPU and x/y-cell of the grid*/
    MPI_Domain_Check(my_size);

    /* (2.2) Domain decomposition */    
    MPI_Dims_create(my_size, BIDIMENSIONAL, dim);
    x_thread = dim[0];
    y_thread = dim[1];


    if (my_rank == 9){
        printf("\n\n");
        printf(" Parallel execution [MPI]\n");
        printf(" =================\n");
        printf("  x-cell: %d \n  y-cell: %d\n", dim[0], dim[1]);
        printf(" =================\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, NO_REORDER, &new_comm);
    MPI_Comm_rank(new_comm, &my_rank);
    
    
    LCPU = MPI_PROC_NULL;
    RCPU = MPI_PROC_NULL;
    TCPU = MPI_PROC_NULL;
    BCPU = MPI_PROC_NULL;
    MPI_Cart_shift(new_comm, 0,  1, &LCPU, &RCPU);
    MPI_Cart_shift(new_comm, 1,  1, &BCPU, &TCPU);

    /* Calcolo la grandezza di ogni cella e ne alloco una per il rispettivo processo ! */
    float *temp_cell;
    int xcell, ycell;
    int *xs, *xe, *ys, *ye;

    

    xs = malloc(my_size*sizeof(int));
    xe = malloc(my_size*sizeof(int));
    ys = malloc(my_size*sizeof(int));
    ye = malloc(my_size*sizeof(int));


    xcell = IMAX/x_thread;
    ycell = JMAX/y_thread;
   

    // total_x_size = ( x_thread == 1 ?    IMAX + 2  :  IMAX + (x_thread - 2)*2 + 2   );
    // total_y_size = ( y_thread == 1 ?    JMAX + 2  :  JMAX + (y_thread - 2)*2 + 2   );

    total_x_size = IMAX + (x_thread)*2 +2 ;
    total_y_size = JMAX + (y_thread)*2 +2 ;

    // size_total_x = size_x+2*x_domains+2;
    // size_total_y = size_y+2*y_domains+2;
    if (my_rank == 0){

        printf("x_thread = %d\n", x_thread);
        printf("y_thread = %d\n", y_thread);
        printf("x_tot = %d\n", total_x_size);
        printf("y_tot = %d\n", total_y_size);
        printf("x_cell = %d\n", xcell);
        printf("y_cell = %d\n", ycell);

    }


    
    /* ================================================================================================= */

    
    processToMap(xs, ys, xe, ye, xcell, ycell, x_thread, y_thread);

    /* Create column data type to communicate with East and West neighBors */
    MPI_Type_vector(xcell,1, JMAX, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);
    
    /* ______________________________________________________________ END */
    
    float **Tn, **Tn1, **Te;
    Tn  = allocateMatrix_MPI(total_x_size, total_y_size);
    Tn1 = allocateMatrix_MPI(total_x_size, total_y_size);
    Te  = allocateMatrix_MPI(IMAX, JMAX);

 
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
        
        initConditionsMPI(Tn,Tn1, total_x_size, total_y_size,x,xD,&time);

        if (my_rank == 0){
            for (i=0; i< total_x_size; i++){
                for (j=0; j< total_y_size; j++){
                    printf("[%d,%d] ", i,j);
                }
                puts("");
            }
        }

        clock_gettime(CLOCK_REALTIME, &time_start);
        
        

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
                if (my_rank == 0){

                    printf("\n\n\t\t\t\tmy_rank = %d\n", my_rank);
                    printf("xe[my_rank] = %d\n", xe[my_rank]);

                    for (i=xs[my_rank];i<=xe[my_rank];i++){
                        for (j=ys[my_rank];j<=ye[my_rank];j++){
                            // printf("%.1f  ", Tn[i][j]);
                            printf("[%d, %d]  ", i-1,j-1);

                        }
                        puts("");
                    }
                }
                MPI_Barrier(new_comm);
              
                if (my_rank == 1){
                    

                    printf("\n\n\t\t\t\tmy_rank = %d\n", my_rank);
                    printf("xe[my_rank] = %d\n", xe[my_rank]);

                    for (i=xs[my_rank];i<=xe[my_rank];i++){
                        for (j=ys[my_rank];j<=ye[my_rank];j++){
                            // printf("%.1f  ", Tn[i][j]);
                            printf("[%d, %d]  ", i-1,j-1);

                        }
                        puts("");
                    }
                }
                MPI_Barrier(new_comm);
                if (my_rank == 2){
                    

                    printf("\n\n\t\t\t\tmy_rank = %d\n", my_rank);
                    printf("xe[my_rank] = %d\n", xe[my_rank]);

                    for (i=xs[my_rank];i<=xe[my_rank];i++){
                        for (j=ys[my_rank];j<=ye[my_rank];j++){
                            printf("[%d, %d]  ", i-1,j-1);

                        }
                        puts("");
                    }
                }
                MPI_Barrier(new_comm);
                if (my_rank == 3){

                    printf("\n\n\t\t\t\tmy_rank = %d\n", my_rank);
                    printf("xe[my_rank] = %d\n", xe[my_rank]);

                    for (i=xs[my_rank];i<=xe[my_rank];i++){
                        for (j=ys[my_rank];j<=ye[my_rank];j++){
                            // printf("%.1f  ", Tn[i][j]);
                            printf("[%d, %d]  ", i-1,j-1);

                        }
                        puts("");
                    }
                }
                MPI_Finalize();
                exit(EXIT_SUCCESS);
                MPI_Barrier(new_comm);
                if (my_rank == 4){

                    printf("\n\n\t\t\t\tmy_rank = %d\n", my_rank);
                    printf("xe[my_rank] = %d\n", xe[my_rank]);

                    for (i=xs[my_rank];i<=xe[my_rank];i++){
                        for (j=ys[my_rank];j<=ye[my_rank];j++){
                            // printf("%.1f  ", Tn[i][j]);
                            printf("[%d, %d]  ", i,j);

                        }
                        puts("");
                    }
                }
                MPI_Barrier(new_comm);
                MPI_Finalize();
                
                exit(EXIT_SUCCESS);
                /* ================================================================================================= */
            
                
                for (i=xs[my_rank];i<=xe[my_rank];i++){
                    for (j=ys[my_rank];j<=ye[my_rank];j++){

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
            



            // updateBound(Tn, RCPU,LCPU,TCPU,BCPU , new_comm, column_type, me, xs, ys, xe, ye, ycell);


            /* Overwrite last solution */
            for (i = xs[my_rank]; i <= xe[my_rank]; i++){
                for (j=ys[my_rank]; j <= ye[my_rank]; j++) {
                    Tn[i][j] = Tn1[i][j];
                }
            }

#ifdef _OUTPUT
            if ( (fabs(time - tio) < 1e-12) ) {

                printf("Plotting data output at time %f\n", time);
                /* Computing exact solution in Te */
                compExactSolution(Te, x, kappa, time);
                DataOutput(test_name, n, time, x, y, Tn, Te);

                plotindex++;
                tio += dtio;
            }
#endif

        }   /* ==== end for (n) ==== */

        clock_gettime(CLOCK_REALTIME, &time_end);

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

