all: int_ring jacobi-mpi
	
int_ring:
	mpicc int_ring.c -o int_ring

jacobi-mpi:
	mpicc jacobi-mpi.c -o jacobi-mpi

clean:
	rm int_ring jacobi-mpi


