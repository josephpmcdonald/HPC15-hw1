#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "util.h"

int main(int argc, char *argv[]) {
    if (2 != argc){
        printf("Incorrect number of args");
        return 0;
    }

    int N = atoi(argv[1]);
    int rank, msgout, msgin, P;
    int tag = 21;
    msgin = 0;
    timestamp_type time1, time2;
    double diff;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    //Get timestamp with processor 0
    int i;
    if (0 == rank) {
        get_timestamp(&time1);
    }

    //Passing single integer as message around the ring
    for (i = 0; i < N; ++i) {

        //Handle processor 0 separately for correct start and end
        if (0 == rank) {
            msgout = msgin + rank;
            //printf("The process %d reporting. The message is %d\n", rank, msgout);
            MPI_Send(&msgout, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
            MPI_Recv(&msgin, 1, MPI_INT, P-1, tag, MPI_COMM_WORLD, &status);
        }
    
        //Other processors
        else {
            MPI_Recv(&msgin, 1, MPI_INT, rank-1, tag, MPI_COMM_WORLD, &status);
            msgout = msgin + rank;
            //printf("The process %d reporting. The message is %d\n", rank, msgout);
            MPI_Send(&msgout, 1, MPI_INT, (rank+1)%P, tag, MPI_COMM_WORLD);
        }
    }

    //Processor 0 gets last timestamp and determines latency
    if (0 == rank) {
        printf("The process %d reporting. The FINAL message is %d\n", rank, msgin);
        get_timestamp(&time2);
        diff = timestamp_diff_in_seconds(time1, time2);
        double latency = diff/(P*N);
        printf("Latency: %f\n", latency);
    }

    //Passing large array around ring
    int size = 250000;
    //double *arrayout = calloc(size, sizeof(double));
    int doublesize = sizeof(double);
    double *array  = calloc(size, doublesize);

    //Processor 0 fills array entries and gets timestamp
    if (0 == rank) {
        for (i = 0; i < size; i++)
            array[i] = i;
         
        get_timestamp(&time1);
    }

    //printf("rank=%d\n",rank);
    for (i = 0; i < N; i++) {
        if (0 == rank) {
            array[0] = i;
            MPI_Send(array, size, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);
            printf("Rank %d, run %f \n", rank, array[0]);
            MPI_Recv(array, size, MPI_DOUBLE, P-1, tag, MPI_COMM_WORLD, &status);
        } 
        else {
            MPI_Recv(array, size, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
            printf("Rank %d, run %f, value in space rank is %f\n", rank, array[0], array[rank]);
            MPI_Send(array, size, MPI_DOUBLE, (rank+1)%P, tag, MPI_COMM_WORLD);
        }
    }

    if (0 == rank) {
        get_timestamp(&time2);
        diff = timestamp_diff_in_seconds(time1, time2);
        //printf("%d %d %d %d %f\n", doublesize, size, P, N, diff);
        //Note that I measure bandwidth as the size of the array in bytes times the number of times it gets passed between processors divided by time. Then reduce it to size in megabytes.
        double bandwidth;
        bandwidth = ((double) doublesize*size*P*N)/(diff*pow(10,6));
        printf("Bandwidth: %f Mbytes per sec\n", bandwidth);
    }

    MPI_Finalize();
    return 0;
}


