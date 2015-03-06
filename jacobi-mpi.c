#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "util.h"


int main(int argc, char* argv[])
{
    int i;
    timestamp_type time1, time2;

    int N = atoi(argv[1]);
    int rank;
    int p;
    MPI_Init();
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (((int) N/p)*p != N) {
        printf("Error: N must be multiple of number of processes\n");
        return 1;
    }

    int n = (int) N/p + 2;

    double *u = (double*) malloc(sizeof(double)*n);
    double *utemp = (double*) malloc(sizeof(double)*n);
    double f = 1.;

    int iter = 0;
    double res = sqrt(N);
    double sum;

    double h =(double) 1/(N+1);
    double hsq = pow(h, 2.);


    for (i = 0; i < N; ++i) {
        u[i] = 0.;
    }

    get_timestamp(&time1);

    while(res > sqrt(N)*1.E-6) {
/*
        if (N > 1000 && iter > 10000)
            break;
*/

        MPI_Send(&u[0], 1, MPI_DOUBLE, rank-1, MPI_COMM_WORLD,

        //Note BC here have u[-1] = 0 for process 0 and u[N] = 0 for last process implicitly
        if (0 == rank):
            utemp[0] = 0;

//            utemp[0] = hsq*(f -(-u[1])/hsq)/2;
        else:
            MPI_Recv
            

        for (i = 1; i < N-1; ++i)
            utemp[i] = hsq*(f - (-u[i-1] - u[i+1])/hsq)/2;
        utemp[N-1] = hsq*(f - (-u[N-2])/hsq)/2;

        for (i = 0; i < N; ++i)
            u[i] = utemp[i];
        
        sum = 0.;
        sum += pow((2*u[0] - u[1])/hsq - 1, 2.);
        for (i = 1; i < N-1; ++i)
            sum += pow((-u[i-1] + 2*u[i] - u[i+1])/hsq - 1, 2.);
        sum += pow((-u[N-2] + 2*u[N-1])/hsq - 1, 2.);
        res = sqrt(sum);
        
        iter++;
    }

    get_timestamp(&time2);
    elapsed = timestamp_diff_in_seconds(time1,time2);

    printf("Jacobi:\n");
    printf("residual = %f\n", res);
    printf("iterations = %d\n", iter);
    printf("elapsed time = %f\n\n", elapsed);
        
    /*
    for (i = 0; i < N; ++i){
        printf("u[%d]=%f\n", i, u[i]);
    }
    */


    free(u);
    free(utemp);

    return 0;
}





