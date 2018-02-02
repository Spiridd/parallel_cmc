#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    // create 3D cartesian grid
    int n_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    int dims[] = {0, 0, 0};
    MPI_Dims_create(n_proc, 3, dims);
    int periods[] = {0, 0, 0};
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &cart_comm);
    // get info about this communicator
    if (cart_comm == MPI_COMM_NULL)
    {
        int oldrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &oldrank);
        if (oldrank == 0)  fprintf(stderr, "cannot create cartesian topology\n");
        exit(EXIT_FAILURE);
    }
    int rank;
    MPI_Comm_rank(cart_comm, &rank);
    int coords[3];
    MPI_Cart_get(cart_comm, 3, dims, periods, coords);
    if (dims[0] != dims[1] || dims[0] != dims[2])
    {
        if (rank == 0)  fprintf(stderr, "n_proc must be cube\n");
        exit(EXIT_FAILURE);
    }
    printf("rank: %d, coord = {%d, %d, %d}\n", rank, coords[0], coords[1], coords[2]);

    // sub communicator
    int dim_subs[] = {0, 0, 1};
    MPI_Comm linear_comm;
    MPI_Cart_sub(cart_comm, dim_subs, &linear_comm);
    if (linear_comm == MPI_COMM_NULL)
    {
        if (rank == 0)  fprintf(stderr, "cannot create linear topology\n");
        exit(EXIT_FAILURE);
    }
    int newrank;
    MPI_Comm_rank(linear_comm, &newrank);
    printf("rank: %d, newrank: %d\n", rank, newrank);

    MPI_Finalize();
    return 0;
}

