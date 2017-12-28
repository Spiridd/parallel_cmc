#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void test_float(int mode, int n, int m, int h, const char* filename_A, const char* filename_B, const char* filename_C)
{
    typedef float value_t;
    #include "fun.def"
    #include "test.def"
}

void test_double(int mode, int n, int m, int h, const char* filename_A, const char* filename_B, const char* filename_C)
{
    typedef double value_t;
    #include "fun.def"
    #include "test.def"
}

void run_float(int n, int m, int h, const char* filename_times)
{
    typedef float value_t;
    #include "fun.def"
    #include "run.def"
}

void run_double(int n, int m, int h, const char* filename_times)
{
    typedef double value_t;
    #include "fun.def"
    #include "run.def"
}

void arg_error(int argc, const char* prog)
{
    printf("Too few arguments.\n" 
            "Usage: %s -(f|d) mode (file_results | file_A file_B file_C) \n"
            "where -f means float, -d means double\n"
            "files are binary and modes are 0, ..., 6.\n"
            "Mode 6 is for report. It requires file_results to be specified.\n", prog);
    exit(1);
}

int main(int argc, char** argv)
{
    if (argc < 4)   arg_error(argc, argv[0]);

    // A = n x m, B = m x h, C = n x h
    static const int n = 500;
    static const int m = 500;
    static const int h = 500;
    
    const char data_t = argv[1][1];
    const int mode = (int) strtol(argv[2], NULL, 10);
    if (mode == 6)
    {
        const char* filename_times = argv[3];
        if (data_t == 'f')  run_float(n, m, h, filename_times);
        else                run_double(n, m, h, filename_times);
    }
    else if (mode<0 || mode>6)
    {
        fprintf(stderr, "incorrect mode\n");
        exit(1);
    }
    else
    {
        if (argc < 6)   arg_error(argc, argv[0]);
        const char* filename_A = argv[3];
        const char* filename_B = argv[4];
        const char* filename_C = argv[5];
        if (data_t == 'f')  test_float(mode, n, m, h, filename_A, filename_B, filename_C);
        else                test_double(mode, n, m, h, filename_A, filename_B, filename_C);
    }

    return 0;
}

