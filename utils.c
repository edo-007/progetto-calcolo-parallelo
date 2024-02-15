#include "utils.h"


char* trim(char *s) {

    char *ptr;
	/* Null string check */
    if ( !s )
        return NULL;   
	/* Empty string check */
    if ( !(*s) )
        return s;      

	/* 
	 * _____isspace()_____________________________________________________
     * If the argument (character) passed to the function is a white-space 
	 * character, it returns non-zero integer. If not, it returns 0.
	*/
    for (ptr = s + strlen(s)-1; (ptr >= s) && isspace(*ptr); --ptr);

    /* Sets string end */
    ptr[1] = '\0';
    return s;
}

float **allocateMatrix() {

    float **matrix = malloc(IMAX * sizeof(float*));

    int i;
    for(i = 0; i < IMAX; i++) {
        matrix[i] = malloc(JMAX*  sizeof(int));
    }

    return matrix;
}

float **allocateMatrix_MPI(int total_x, int total_y) {

    float **matrix = malloc(total_x * sizeof(float*));

    int i;
    for(i = 0; i < total_x; i++) {
        matrix[i] = malloc(total_y*  sizeof(int));
    }

    return matrix;
}

void initConditions( float **Tn,float **Tn1, float **Te, float *x, float xD,float *time){
    int i,j;
    for (i = 0; i < IMAX; i++){
        for (j = 0; j < JMAX; j++){
                // printf("(%d, %d) ", i,j);
                if ( x[i] < xD )
                    Tn[i][j] = TL;
                else
                    Tn[i][j] = TR;    
                
                Tn1[i][j] = 0.0;
        }
        // puts("");
    }
    *time = (float)0;
}

void compExactSolution(float **Te, float *x, float kappa, float time){
    
    int i,j; 
    for (i = 0; i < IMAX; i++) {
        for(j = 0; j < JMAX; j++) {
            Te[i][j] = Te[i][j] = ((TR + TL) / 2.0) + (  erf( (x[i] - 1) / (2 * sqrt(kappa * time)) ) * (TR - TL) / 2.0  );
        }
    }
}

void DataOutput(  char testname[200], int timestep, float time, float* x, float* y, float **T, float **Te) {

    int i,j, DataUnit;
    char citer[10];
    char IOFilename[200];
	char error_mex[250];
    FILE *fp;

    trim(testname);
    sprintf(IOFilename, "%s/%s-%04d.csv", DAT_PATH, testname, timestep);

    if ( ( fp = fopen( IOFilename, "w+t") ) == NULL ) {
        /* Errore nell' apertura del file */
        sprintf(error_mex, "errore nell' apertura del file %s", IOFilename ); 
        perror(error_mex);
        exit(EXIT_FAILURE);
    }

    /* Scrittura su file */
    fprintf(fp, "CURRENT TIME : %f\n", time );
    fprintf(fp, "VARIABLES : y - x - T - Te\n" );
    fprintf(fp, "ZONE T='Only Zone', I= %d, J=%d  F=POINT\n", IMAX, JMAX  );
    
    for (i = 0; i < IMAX; i++) {
        for (j = 0; j < JMAX; j++) {
            fprintf(fp, "%.15f, %.15f, %.15f, %.15f\n", y[i], x[j], T[i][j], Te[i][j]);
        }
    }

    fclose(fp);

 }

void Time_Output( float elapsed_time ) {

    FILE* fp;

    if ( ( fp = fopen("times.csv", "w") ) == NULL)
    {
        perror("can' t open time-output file");
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "%.8f, 1", elapsed_time);

    fclose(fp);
}


void sub_timespec ( struct timespec t1, struct timespec t2, struct timespec *td) {

	td->tv_nsec = t2.tv_nsec - t1.tv_nsec;
	td->tv_sec = t2.tv_sec - t1.tv_sec;

	
	if (td->tv_sec > 0 && td->tv_nsec < 0 ){
		td->tv_nsec += 	NS_PER_SECOND;
		td->tv_sec--;
	}
	if (td->tv_sec < 0 && td->tv_nsec > 0 ){
                td->tv_nsec -=  NS_PER_SECOND;
                td->tv_sec++;
        }
}

double simple_sub_timespec ( struct timespec t1 , struct timespec t2 ) {

	double td1, td2 ;

	td1 = t1.tv_sec + ( t1.tv_nsec /( double ) NS_PER_SECOND ) ;
	td2 = t2.tv_sec + ( t2.tv_nsec /( double ) NS_PER_SECOND ) ;

	return ( td2 - td1 ) ;
}

void print_tn(float **mat){

    int i,j;
    for (i=0;i< IMAX; i++){
        for (j=0; j< JMAX; j++){
            printf("%.1f ", mat[i][j] );
        }
        printf("\n");
    }
    
}

