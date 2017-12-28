#include <stdio.h>
#include <stdlib.h>
#include <papi.h>
#include "matrix.h"

void error_handler(const char* str)
{
    fprintf(stderr, "%s", str);
    exit(1);
}

// measures and prints to file times and PAPI events counts
void test(unsigned int* Events, int NUM_EVENTS, FILE* file)
{
    const static int nmin = 1000;
    const static int nmax = 5000;
    const static int step = 1000;

    float* A = malloc(nmax*nmax*sizeof(float));
    float* B = malloc(nmax*nmax*sizeof(float));
    float* C = malloc(nmax*nmax*sizeof(float));
    if (A==NULL || B==NULL || C==NULL)
        error_handler("cannot malloc");

    int EventSet = PAPI_NULL;
    if (PAPI_create_eventset(&EventSet) != PAPI_OK)
        error_handler("cannot create eventset\n");
    if (PAPI_add_events(EventSet, Events, NUM_EVENTS) != PAPI_OK)
        error_handler("cannot add events\n");
    if (PAPI_start(EventSet) != PAPI_OK)
        error_handler("cannot start papi\n");
    long long values[NUM_EVENTS];

    // all sizes
    for(int n=nmin; n<=nmax; n+=step)
    {
        printf("size = %d\n", n);
        // two runs: 32x32 and 52x52
        for(int bsz=32; bsz<=52; bsz+=20)
        {
            // ijk ijk
            set_zeros(C, n*n);
            get_random_array(A, n*n);
            get_random_array(B, n*n);

            long long begin = PAPI_get_virt_usec();
            if (PAPI_reset(EventSet) != PAPI_OK)
                error_handler("cannot reset counters\n");
            square_dgemm_ijk_with_call_ijk(A, B, C, n, bsz);
            if (PAPI_read(EventSet, values) != PAPI_OK)
                error_handler("cannot read values\n");
            long long end = PAPI_get_virt_usec();

            fprintf(file, "ijk,%d,%d,", n, bsz);
            for(int i=0; i<NUM_EVENTS; ++i)
                fprintf(file, "%llu,", values[i]);
            double time = (double) (end-begin)/1e6;
            fprintf(file, "%f\n", time);

            // ijk ikj - 32 only
            if (bsz != 32)  continue;
            set_zeros(C, n*n);
            get_random_array(A, n*n);
            get_random_array(B, n*n);

            begin = PAPI_get_virt_usec();
            if (PAPI_reset(EventSet) != PAPI_OK)
                error_handler("cannot reset counters\n");
            square_dgemm_ijk_with_call_ikj(A, B, C, n, bsz);
            if (PAPI_read(EventSet, values) != PAPI_OK)
                error_handler("cannot read values\n");
            end = PAPI_get_virt_usec();

            fprintf(file, "ikj,%d,%d,", n, bsz);
            time = (double) (end-begin)/1e6;
            for(int i=0; i<NUM_EVENTS; ++i)
                fprintf(file, "%llu,", values[i]);
            fprintf(file, "%f\n", time);
        }
    }

    if (PAPI_stop(EventSet, values) != PAPI_OK)
        error_handler("cannot stop papi\n");

    free(A);
    free(B);
    free(C);
}

int main()
{
    // PAPI initialization
    if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
        error_handler("cannot init\n");

    FILE* file1 = fopen("../res/results1.txt", "w");
    fprintf(file1, "indices,size,bsz,PAPI_L1_DCM,PAPI_L2_DCM,PAPI_TOT_CYC,time\n");
    unsigned int Events1[] = {PAPI_L1_DCM, PAPI_L2_DCM, PAPI_TOT_CYC};
    test(Events1, sizeof(Events1)/sizeof(Events1[0]), file1);
    fclose(file1);

    // Events2 can not be measured simultaneously with Events1
    FILE* file2 = fopen("../res/results2.txt", "w");
    fprintf(file2, "indices,size,bsz,PAPI_FP_OPS,PAPI_TLB_IM,PAPI_TLB_DM,time\n");
    unsigned int Events2[] = {PAPI_FP_OPS, PAPI_TLB_IM, PAPI_TLB_DM};
    test(Events1, sizeof(Events2)/sizeof(Events2[0]), file2);
    fclose(file2);


    return 0;
}

