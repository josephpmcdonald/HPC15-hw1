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
    int maxiter = atoi(argv[2]);

    int rank;
    int p;
    int tag = 21;
    MPI_Init();
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (((int) N/p)*p != N) {
        printf("Error: N must be multiple of number of processes\n");
        return 1;
    }

    int n = (int) N/p + 2;

    double h =(double) 1/(N+1);
    double hsq = pow(h, 2.);

    double *u = (double*) malloc(sizeof(double)*n);
    double *utemp = (double*) malloc(sizeof(double)*n);
    double f = 1.;

    int iter = 0;
    double res = sqrt(N);
    double sum;
    double other;//To collect residuals from processes other than 0


    for (i = 0; i < n; ++i) {
        u[i] = 0.;
    }

    get_timestamp(&time1);

    while(res > sqrt(N)*1.E-6 && iter < maxiter) {
        //BC here have u[0] = 0 for process 0
        if (0 == rank) {
            u[0] = 0;
            MPI_Send(&u[n-2], 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
            MPI_Recv(&u[n-1], 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        }

        else if (p-1 == rank) {
            u[n-1] = 0;
            MPI_Recv(&u[0], 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
            MPI_Send(&u[1], 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
        }

        else {
            MPI_Recv(&u[0], 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
            MPI_Send(&u[1], 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
            MPI_Send(&u[n-2], 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
            MPI_Recv(&u[n-1], 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status); 
        }


//        utemp[0] = hsq*(f -(-u[1])/hsq)/2;

        for (i = 1; i < n-1; ++i)
            utemp[i] = hsq*(f - (-u[i-1] - u[i+1])/hsq)/2;

//        utemp[N-1] = hsq*(f - (-u[N-2])/hsq)/2;

        for (i = 1; i < n-1; ++i)
            u[i] = utemp[i];
        
        sum = 0.;

        /*Note we don't need to add to the residual for the first cell in process 0
          or the last cell in process p-1 since they're set by the BC.*/
        for (i = 1; i < n-1; ++i)
            sum += pow((-u[i-1] + 2*u[i] - u[i+1])/hsq - f, 2.);

        /*Receive the residual sums from all other processes, add them together
          take square root and send back to other processes*/
        if (0 == rank) {
            for (i = 1; i < p; ++i) {
                MPI_Recv(&other, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
                sum += other;
            }
            
            res = sqrt(sum);
            for (i = 1; i < p; ++i)
                MPI_Send(&res, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
        }
        //Receive the residual from process 0
        else {
            MPI_Send(&sum, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
            MPI_Recv(&res, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        }

        
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


//    free(u);
//    free(utemp);

    return 0;
}





