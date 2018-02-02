#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include "fun.h"

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    if (argc != 4)
    {
        fprintf(stderr, "Incorrect run.\n"
                "Usage: %s file_A file_b file_c\n"
                "where file_A is name of file"
                " where matrix A is located.\n", argv[0]);
        exit(1);
    }

    parallel_multiply(argv[1], argv[2], argv[3]);

    MPI_Finalize();

    return 0;
}

