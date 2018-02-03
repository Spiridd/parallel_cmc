/*
 * Generetes random matrix and random vector,
 * puts them in file1 and file2
 */
#include <stdio.h>
#include <stdlib.h>
#include "fun.h"

int main(int argc, char** argv)
{
    if (argc != 5)
    {
        fprintf(stderr, "Incorrect args.\n"
                "Usage: %s m n file1 file2\n"
                "where m is #rows, n is #coloumns;\n"
                "file1 is for matrix, file2 is for vector.\n", argv[0]);
        exit(1);
    }
    const int m = (int) strtol(argv[1], NULL, 10);
    const int n = (int) strtol(argv[2], NULL, 10);
    if (n<0 || m<0)
    {
        fprintf(stderr, "Provide positive sizes\n");
        exit(1);
    }
    if (n != m)
    {
        fprintf(stderr, "N must be equal to M\n");
        exit(EXIT_FAILURE);
    }

    double* const M = safe_malloc(n*n*sizeof(double));
    get_random_array(M, n*n);
    write_array_to_file(M, n, n, argv[3]);

    double* const V = safe_malloc(n*n*sizeof(double));
    get_random_array(V, n*n);
    write_array_to_file(V, n, n, argv[4]);

    return 0;
}

