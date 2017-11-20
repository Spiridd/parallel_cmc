#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix.h"

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
    int mode = (int) strtol(argv[4], NULL, 10);
    char *filename_times = NULL;
    int is_test = 1;
    if (mode == 0)
    {
        filename_times = argv[5];
        is_test = 0;
    }

    const int n = 500;
    double **A = allocate_matrix_2D(n);
    double **B = allocate_matrix_2D(n);
    double **C = allocate_matrix_2D(n);

    get_random_matrix_2D(A, n);
    get_random_matrix_2D(B, n);
    set_zeros_2D(C, n);

    void (*fp[6])(double**, double**, double**, int);

    fp[0] = &multiply_ijk_2D;
    fp[1] = &multiply_jik_2D;
    fp[2] = &multiply_ikj_2D;
    fp[3] = &multiply_kij_2D;
    fp[4] = &multiply_kji_2D;
    fp[5] = &multiply_jki_2D;

    if (is_test)
    {
        if (mode <= 6)
        {
            fp[mode-1](A, B, C, n);
            write_matrix_to_binary_2D(A, n, filename_A);
            write_matrix_to_binary_2D(B, n, filename_B);
            write_matrix_to_binary_2D(C, n, filename_C);
        }
        else
        {
            printf("Incorrect run");
            return 1;
        }
    }
    else
    {
        const int n_iter = 10;
        double times[6];

        for(int i=0; i<6; ++i)
            times[i] = measure_time_2D(fp[i], A, B, C, n, n_iter);
        
        write_time_to_binary_2D(times, sizeof(times)/sizeof(times[0]), filename_times);
    }

    free_matrix_2D(A, n);
    free_matrix_2D(B, n);
    free_matrix_2D(C, n);

    return 0;
}

