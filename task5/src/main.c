#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "matrix.h"

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    if (argc != 4)
    {
        fprintf(stderr, "Incorrect run.\n"
                "Usage: %s file_A file_B file_C\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    parallel_multiply(argv[1], argv[2], argv[3]);

    MPI_Finalize();

    return 0;
}

