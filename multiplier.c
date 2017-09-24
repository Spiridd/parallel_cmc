#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void write_time_to_binary(double times[], size_t count, char *filename)
{
    FILE* file = fopen(filename, "wb");
    fwrite(times, sizeof(double), count, file);
    fclose(file);
}

void write_matrix_to_binary(double **A, int n, char *filename)
{
    FILE *file = fopen(filename, "wb");
    for(int i=0; i<n; ++i)
        fwrite(A[i], sizeof(double), n, file);
    fclose(file);
}

void print_matrix(double **A, int n)
{
    for(int i=0; i<n; ++i)
    {
        for(int j=0; j<n; ++j)
            printf("%lf ", A[i][j]);
        printf("\n");
    }
}

void get_random_matrix(double **M, int n)
{
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            M[i][j] = (double)rand()/RAND_MAX;
}

void set_zeros(double **A, int n)
{
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            A[i][j] = 0.0;
}

double** allocate_matrix(int n)
{
    double **A = malloc(n * sizeof(*A));
    for(int i=0; i<n; ++i)
        A[i] = malloc(n * sizeof(*(A[i])));
    return A;
}

void free_matrix(double **A, int n)
{
    for(int i=0; i<n; ++i)
        free(A[i]);
    free(A);
}

void multiply_1(double **A, double **B, double **C, int n)
{
    // matrix C is filled in row order
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
        {
            for(int k=0; k<n; ++k)
                C[i][j] += A[i][k] * B[k][j];
        }
}

void multiply_2(double **A, double **B, double **C, int n)
{
    // matrix C is filled in column order
    for(int j=0; j<n; ++j)
        for(int i=0; i<n; ++i)
        {
            for(int k=0; k<n; ++k)
                C[i][j] += A[i][k] * B[k][j];
        }
}
void multiply_3(double **A, double **B, double **C, int n)
{
    // matrix A is used in row order
    for(int i=0; i<n; ++i)
        for(int k=0; k<n; ++k)
            for(int j=0; j<n; ++j)
                C[i][j] += A[i][k] * B[k][j];
}
void multiply_4(double **A, double **B, double **C, int n)
{
    // matrix A is used in column order
    for(int k=0; k<n; ++k)
        for(int i=0; i<n; ++i)
            for(int j=0; j<n; ++j)
                C[i][j] += A[i][k] * B[k][j];
}
void multiply_5(double **A, double **B, double **C, int n)
{
    // matrix B is used in row order
    for(int k=0; k<n; ++k)
        for(int j=0; j<n; ++j)
            for(int i=0; i<n; ++i)
                C[i][j] += A[i][k] * B[k][j];
}
void multiply_6(double **A, double **B, double **C, int n)
{
    // matrix B is used in column order
    for(int j=0; j<n; ++j)
        for(int k=0; k<n; ++k)
            for(int i=0; i<n; ++i)
                C[i][j] += A[i][k] * B[k][j];
}

int main(int argc, char** argv)
{
    if (argc < 5)
    {
        printf("Too few arguments." 
                "Usage: %s file_A file_B file_C mode [file_results]\n"
                "where files are binary and mode is 0, ..., 6.\n"
                "Mode 0 requires file_results to be specified.\n", argv[0]);
        return 1;
    }

    srand(time(NULL));
    
    char *filename_A = argv[1];
    char *filename_B = argv[2];
    char *filename_C = argv[3];
    long mode = strtol(argv[4], NULL, 10);
    char *filename_times = NULL;
    int is_test = 1;
    if (mode == 0)
    {
        filename_times = argv[5];
        is_test = 0;
    }

    const int n = 300;
    double **A = allocate_matrix(n);
    double **B = allocate_matrix(n);
    double **C = allocate_matrix(n);

    get_random_matrix(A, n);
    //print_matrix(A, n);
    get_random_matrix(B, n);
    set_zeros(C, n);

    if (is_test)
    {
        switch (mode)
        {
            case 1:
                multiply_1(A, B, C, n);
                break;
            case 2:
                multiply_1(A, B, C, n);
                break;
            case 3:
                multiply_1(A, B, C, n);
                break;
            case 4:
                multiply_1(A, B, C, n);
                break;
            case 5:
                multiply_1(A, B, C, n);
                break;
            case 6:
                multiply_1(A, B, C, n);
                break;
            default:
                printf("Incorrect run");
                return 1;
        }
        write_matrix_to_binary(A, n, filename_A);
        write_matrix_to_binary(B, n, filename_B);
        write_matrix_to_binary(C, n, filename_C);
    }
    else
    {
        const int n_iter = 10;
        double times[6]; // six functions
        clock_t begin, end;

        begin = clock();
        for(int i=0; i<n_iter; ++i)
            multiply_1(A, B, C, n);
        end = clock();
        times[0] = (double) (end-begin)/CLOCKS_PER_SEC/n_iter;

        begin = clock();
        for(int i=0; i<n_iter; ++i)
            multiply_2(A, B, C, n);
        end = clock();
        times[1] = (double) (end-begin)/CLOCKS_PER_SEC/n_iter;

        begin = clock();
        for(int i=0; i<n_iter; ++i)
            multiply_3(A, B, C, n);
        end = clock();
        times[2] = (double) (end-begin)/CLOCKS_PER_SEC/n_iter;

        begin = clock();
        for(int i=0; i<n_iter; ++i)
            multiply_4(A, B, C, n);
        end = clock();
        times[3] = (double) (end-begin)/CLOCKS_PER_SEC/n_iter;

        begin = clock();
        for(int i=0; i<n_iter; ++i)
            multiply_5(A, B, C, n);
        end = clock();
        times[4] = (double) (end-begin)/CLOCKS_PER_SEC/n_iter;

        begin = clock();
        for(int i=0; i<n_iter; ++i)
            multiply_6(A, B, C, n);
        end = clock();
        times[5] = (double) (end-begin)/CLOCKS_PER_SEC/n_iter;
        
        write_time_to_binary(times, sizeof(times)/sizeof(times[0]), filename_times);
    }

    free_matrix(A, n);
    free_matrix(B, n);
    free_matrix(C, n);

    return 0;
}
