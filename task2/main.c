#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "matrix.h"

void test_block_and_ikj(int n, int bsz)
{
    double *A = allocate_matrix_1D(n);
    double *B = allocate_matrix_1D(n);
    double *C1 = allocate_matrix_1D(n);
    double *C2 = allocate_matrix_1D(n);

    get_random_matrix_1D(A, n);
    get_random_matrix_1D(B, n);
    set_zeros_1D(C1, n);
    set_zeros_1D(C2, n);

    void (*f1)(double*, double*, double*, int) = &multiply_ikj_1D;
    void (*f2)(double*, double*, double*, int, int) = &multiply_block_1D;
    //void (*f2)(double*, double*, double*, int, int) = &square_dgemm_1D;

    f1(A, B, C1, n);
    f2(A, B, C2, n, bsz);

    static const double eps = 1e-6;
    int errors = 0;
    for(int i=0; i<n*n; ++i)
        if (abs(C1[i]-C2[i])>eps)       ++errors;
    printf("there are %d positions and %d errors\n", n*n, errors);
        
    free_matrix_1D(A, n);
    free_matrix_1D(B, n);
    free_matrix_1D(C1, n);
    free_matrix_1D(C2, n);

}

void time_all_2D(int n, int n_iter)
{
    double **A = allocate_matrix_2D(n);
    double **B = allocate_matrix_2D(n);
    double **C = allocate_matrix_2D(n);

    get_random_matrix_2D(A, n);
    get_random_matrix_2D(B, n);
    set_zeros_2D(C, n);

    static const int n_func = 6;
    void (*fp[n_func])(double**, double**, double**, int);
    fp[0] = &multiply_ijk_2D;
    fp[1] = &multiply_jik_2D;
    fp[2] = &multiply_ikj_2D;
    fp[3] = &multiply_kij_2D;
    fp[4] = &multiply_kji_2D;
    fp[5] = &multiply_jki_2D;

    double times[n_func];
    for(int i=0; i<n_func; ++i)
    {
        times[i] = measure_time_2D(fp[i], A, B, C, n, n_iter);
        printf("%f\n", times[i]);
    }
        
    free_matrix_2D(A, n);
    free_matrix_2D(B, n);
    free_matrix_2D(C, n);
}

void time_all_1D(int n, int n_iter)
{
    double *A = allocate_matrix_1D(n);
    double *B = allocate_matrix_1D(n);
    double *C = allocate_matrix_1D(n);

    get_random_matrix_1D(A, n);
    get_random_matrix_1D(B, n);
    set_zeros_1D(C, n);

    static const int n_func = 6;
    void (*fp[n_func])(double*, double*, double*, int);
    fp[0] = &multiply_ijk_1D;
    fp[1] = &multiply_jik_1D;
    fp[2] = &multiply_ikj_1D;
    fp[3] = &multiply_kij_1D;
    fp[4] = &multiply_kji_1D;
    fp[5] = &multiply_jki_1D;

    double times[n_func];
    for(int i=0; i<n_func; ++i)
    {
        times[i] = measure_time_1D(fp[i], A, B, C, n, n_iter);
        printf("%f\n", times[i]);
    }
        
    free_matrix_1D(A, n);
    free_matrix_1D(B, n);
    free_matrix_1D(C, n);
}

void time_block_1D(int n, int n_iter, int bsz)
{
    double *A = allocate_matrix_1D(n);
    double *B = allocate_matrix_1D(n);
    double *C = allocate_matrix_1D(n);

    get_random_matrix_1D(A, n);
    get_random_matrix_1D(B, n);
    set_zeros_1D(C, n);

    //void (*fp)(double*, double*, double*, int, int) = &multiply_block_1D;
    void (*fp)(double*, double*, double*, int, int) = &square_dgemm_1D;
    clock_t begin, end;
    begin = clock();
    for(int i=0; i<n_iter; ++i)
        fp(A, B, C, n, bsz);
    end = clock();
    const double time = (double)(end-begin)/CLOCKS_PER_SEC/n_iter;
    printf("block size = %d, time = %f s\n", bsz, time);
        
    free_matrix_1D(A, n);
    free_matrix_1D(B, n);
    free_matrix_1D(C, n);
}

int main(int argc, char** argv)
{
    static const int n = 1000; // matrix size
    printf("matrix size = %dx%d\n", n, n);
    static const int n_iter = 3;
    //static const int bsz = 36;

    //printf("matrix as 2D array:\n");
    //time_all_2D(n, n_iter);

    //printf("matrix as 1D array:\n");
    //time_all_1D(n, n_iter);
    
    printf("block algorithm:\n");
    for(int bsz=10; bsz<50; ++bsz)
    {
        time_block_1D(n, n_iter, bsz);
        test_block_and_ikj(n, bsz);
    }

    return 0;
}

