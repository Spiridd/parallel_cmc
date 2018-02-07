#include <stdlib.h>
#include <stdio.h>
#include "fun.h"

void write_array_to_file(double *A, int m, int n, const char *filename)
{
    FILE *file = fopen(filename, "wb");
    // 0 - float, 1 - double
    // here we use double only
    // (write policy is compatible with task1)
    const static int f_type = 0;
    const static int d_type = 1;
    fwrite(&d_type, sizeof(int), 1, file);
    fwrite(&m, sizeof(int), 1, file);
    fwrite(&n, sizeof(int), 1, file);
    fwrite(A, sizeof(double), m*n, file);
    fclose(file);
}

void get_random_array(double *A, int N)
{
    for(int i=0; i<N; ++i)
        A[i] = (double)rand()/RAND_MAX;
}

void* safe_malloc(size_t nbytes)
{
    void *ret = malloc(nbytes);
    if (ret==NULL)
    {
        fprintf(stderr, "cannot malloc\n");
        exit(1);
    }
    return ret;
}

