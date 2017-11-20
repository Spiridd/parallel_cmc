#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix.h"

void write_time_to_binary_2D(double times[], size_t count, char *filename)
{
    FILE* file = fopen(filename, "wb");
    fwrite(times, sizeof(double), count, file);
    fclose(file);
}

void write_matrix_to_binary_2D(double **A, int n, char *filename)
{
    FILE *file = fopen(filename, "wb");
    for(int i=0; i<n; ++i)
        fwrite(A[i], sizeof(double), n, file);
    fclose(file);
}

void print_matrix_2D(double **A, int n)
{
    for(int i=0; i<n; ++i)
    {
        for(int j=0; j<n; ++j)
            printf("%lf ", A[i][j]);
        printf("\n");
    }
}

void get_random_matrix_2D(double **M, int n)
{
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            M[i][j] = (double)rand()/RAND_MAX;
}

void set_zeros_2D(double **A, int n)
{
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            A[i][j] = 0.0;
}

double** allocate_matrix_2D(int n)
{
    double **A = malloc(n * sizeof(*A));
    for(int i=0; i<n; ++i)
        A[i] = malloc(n * sizeof(*(A[i])));

    return A;
}

void free_matrix_2D(double **A, int n)
{
    for(int i=0; i<n; ++i)
        free(A[i]);
    
    free(A);
}

void multiply_ijk_2D(double **A, double **B, double **C, int n)
{
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
        {
            double sum = 0.0;
            for(int k=0; k<n; ++k)
                sum += A[i][k] * B[k][j];
            C[i][j] = sum;
        }
}

void multiply_jik_2D(double **A, double **B, double **C, int n)
{
    for(int j=0; j<n; ++j)
        for(int i=0; i<n; ++i)
        {
            double sum = 0.0;
            for(int k=0; k<n; ++k)
                sum += A[i][k] * B[k][j];
            C[i][j] = sum;
        }
}

void multiply_ikj_2D(double **A, double **B, double **C, int n)
{
    for(int i=0; i<n; ++i)
        for(int k=0; k<n; ++k)
        {
            double a = A[i][k];
            for(int j=0; j<n; ++j)
                C[i][j] += a * B[k][j];
        }
}

void multiply_kij_2D(double **A, double **B, double **C, int n)
{
    for(int k=0; k<n; ++k)
        for(int i=0; i<n; ++i)
        {
            double a = A[i][k];
            for(int j=0; j<n; ++j)
                C[i][j] += a * B[k][j];
        }
}

void multiply_kji_2D(double **A, double **B, double **C, int n)
{
    for(int k=0; k<n; ++k)
        for(int j=0; j<n; ++j)
            for(int i=0; i<n; ++i)
            {
                double b = B[k][j];
                C[i][j] += A[i][k] * B[k][j];
            }
}

void multiply_jki_2D(double **A, double **B, double **C, int n)
{
    for(int j=0; j<n; ++j)
        for(int k=0; k<n; ++k)
            for(int i=0; i<n; ++i)
            {
                double b = B[k][j];
                C[i][j] += A[i][k] * B[k][j];
            }
}

double measure_time_2D(void (*multiplier)(double**, double**, double**, int), double **A, double **B, double **C, int n, int n_iter)
{
    clock_t begin, end;

    begin = clock();
    for(int i=0; i<n_iter; ++i)
        multiplier(A, B, C, n);
    end = clock();

    return (double)(end-begin)/CLOCKS_PER_SEC/n_iter;
}
