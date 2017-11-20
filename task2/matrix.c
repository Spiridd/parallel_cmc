#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix.h"
#define min(a, b) (a < b ? a : b)

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

void get_random_matrix_1D(double *M, int n)
{
    for(int i=0; i<n; ++i)
    {
        int i_offset = i*n;
        for(int j=0; j<n; ++j)
            M[i_offset + j] = (double)rand()/RAND_MAX;
    }
}

void set_zeros_2D(double **A, int n)
{
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            A[i][j] = 0.0;
}

void set_zeros_1D(double *A, int n)
{
    for(int i=0; i<n; ++i)
    {
        int i_offset = i*n;
        for(int j=0; j<n; ++j)
            A[i_offset + j] = 0.0;
    }
}

double** allocate_matrix_2D(int n)
{
    double **A = malloc(n * sizeof(*A));
    for(int i=0; i<n; ++i)
        A[i] = malloc(n * sizeof(*(A[i])));

    return A;
}

double* allocate_matrix_1D(int n)
{
    return malloc(n*n * sizeof(double));
}

void free_matrix_2D(double **A, int n)
{
    for(int i=0; i<n; ++i)
        free(A[i]);
    
    free(A);
}

void free_matrix_1D(double *A, int n)
{
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
        {
            const double b = B[k][j];
            for(int i=0; i<n; ++i)
                C[i][j] += A[i][k] * b;
        }
}

void multiply_jki_2D(double **A, double **B, double **C, int n)
{
    for(int j=0; j<n; ++j)
        for(int k=0; k<n; ++k)
        {
            const double b = B[k][j];
            for(int i=0; i<n; ++i)
                C[i][j] += A[i][k] * b;
        }
}

void multiply_ijk_1D(double *A, double *B, double *C, int n)
{
    for(int i=0; i<n; ++i)
    {
        const int i_offset = i*n;
        for(int j=0; j<n; ++j)
        {
            double sum = C[i_offset + j];
            for(int k=0; k<n; ++k)
                sum += A[i_offset + k] * B[k*n + j];
            C[i_offset + j] = sum;
        }
    }
}

void multiply_jik_1D(double *A, double *B, double *C, int n)
{
    for(int j=0; j<n; ++j)
        for(int i=0; i<n; ++i)
        {
            const int i_offset;
            double sum = C[i_offset + j];
            for(int k=0; k<n; ++k)
                sum += A[i_offset + k] * B[k*n + j];
            C[i_offset + j] = sum;
        }
}

void multiply_ikj_1D(double *A, double *B, double *C, int n)
{
    for(int i=0; i<n; ++i)
    {
        const int i_offset = i*n;
        for(int k=0; k<n; ++k)
        {
            const int k_offset = k*n;
            double a = A[i_offset + k];
            for(int j=0; j<n; ++j)
                C[i_offset + j] += a * B[k_offset + j];
        }
    }
}

void multiply_kij_1D(double *A, double *B, double *C, int n)
{
    for(int k=0; k<n; ++k)
    {
        const int k_offset = k*n;
        for(int i=0; i<n; ++i)
        {
            const int i_offset = i*n;
            double a = A[i_offset + k];
            for(int j=0; j<n; ++j)
                C[i_offset + j] += a * B[k_offset + j];
        }
    }
}

void multiply_kji_1D(double *A, double *B, double *C, int n)
{
    for(int k=0; k<n; ++k)
    {
        const int k_offset = k*n;
        for(int j=0; j<n; ++j)
        {
            const double b = B[k_offset + j];
            for(int i=0; i<n; ++i)
                C[i*n + j] += A[i*n + k] * b;
        }
    }
}

void multiply_jki_1D(double *A, double *B, double *C, int n)
{
    for(int j=0; j<n; ++j)
        for(int k=0; k<n; ++k)
        {
            const int k_offset = k*n;
            const double b = B[k_offset + j];
            for(int i=0; i<n; ++i)
                C[i*n + j] += A[i*n + k] * b;
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

double measure_time_1D(void (*multiplier)(double*, double*, double*, int), double *A, double *B, double *C, int n, int n_iter)
{
    clock_t begin, end;

    begin = clock();
    for(int i=0; i<n_iter; ++i)
        multiplier(A, B, C, n);
    end = clock();

    return (double)(end-begin)/CLOCKS_PER_SEC/n_iter;
}

void calculate_block_1D(double *A, double *B, double *C, int M, int N, int K, int n)
{
    for(int i=0; i<M; ++i)
    {
        const int i_offset = i*n;
        for(int j=0; j<N; ++j)
        {
            double cij = 0.0;
            int k=0; 
            for( ; k<K-8; k+=8)
            {
                const double d0 = A[i_offset+k] * B[k*n+j];
                const double d1 = A[i_offset+k+1] * B[(k+1)*n+j];
                const double d2 = A[i_offset+k+2] * B[(k+2)*n+j];
                const double d3 = A[i_offset+k+3] * B[(k+3)*n+j];
                const double d4 = A[i_offset+k+4] * B[(k+4)*n+j];
                const double d5 = A[i_offset+k+5] * B[(k+5)*n+j];
                const double d6 = A[i_offset+k+6] * B[(k+6)*n+j];
                const double d7 = A[i_offset+k+7] * B[(k+7)*n+j];
                cij += (d0+d1+d2+d3+d4+d5+d6+d7);
            }
            
            for( ; k<K; ++k)
                cij += A[i_offset+k] * B[k*n+j];

            C[i_offset+j] += cij;
        }
    }
}

void multiply_block_1D(double *A, double *B, double *C, int n, int bsz)
{
    for(int i=0; i<n; i+=bsz)
    {
        const int i_offset = i*n;
        for(int j=0; j<n; j+=bsz)
            for(int k=0; k<n; k+=bsz)
            {
                // correcting case of the last block
                const int M = min(bsz, n-i);
                const int N = min(bsz, n-j);
                const int K = min(bsz, n-k);
                // block multiplication
                for(int i1=0; i1<M; ++i1)
                {
                    const int i1_offset = i1*n;
                    for(int j1=0; j1<N; ++j1)
                    {
                        double cij = 0.0;
                        for(int k1=0; k1<K; ++k1)
                            cij += A[i_offset+i1_offset+k+k1]
                                   * B[(k+k1)*n+j+j1];
                        C[i_offset+i1_offset+j+j1] += cij;
                    }
                }
            }
    }
}

void square_dgemm_1D(double *A, double *B, double *C, int n, int bsz)
{
    for(int i=0; i<n; i+=bsz)
    {
        const int i_offset = i*n;
        for(int j=0; j<n; j+=bsz)
            for(int k=0; k<n; k+=bsz)
            {
                const int M = min(bsz, n-i);
                const int N = min(bsz, n-j);
                const int K = min(bsz, n-k);

                calculate_block_1D(A+i_offset+k, B+k*n+j, C+i_offset+j, M, N, K, n);
            }
    }
}

