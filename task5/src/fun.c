#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
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

void distribute_tasks(int n_tasks, int n_proc, int **tasks_ptr, int **from_ptr, int **to_ptr)
{
    /* distributing tasks in a following manner:
     * (4 processes, 22 tasks)
     * 0    1   2   3   ranks
     *
     * 5    5   5   5   not all tasks
     * 6    6   5   5   all the tasks (with remained)
     */
    const int tpp = (n_tasks-1)/n_proc + 1; // max tasks per process
    const int remained = n_tasks - (tpp-1)*n_proc;
    const int min_task_quantity = (tpp > 0 ? tpp-1 : 0);
    *tasks_ptr = safe_malloc(n_proc*sizeof(int));
    int * const tasks = *tasks_ptr;
    for(int i=0; i<n_proc; ++i) tasks[i] = min_task_quantity;
    for(int i=0; i<remained; ++i) ++tasks[i];

    *from_ptr = safe_malloc(n_proc*sizeof(int));
    *to_ptr = safe_malloc(n_proc*sizeof(int));
    int * const from = *from_ptr;
    int * const to = *to_ptr;
    from[0] = 0;
    to[0] = tasks[0];
    for(int i=1; i<n_proc; ++i)
    {
        from[i] = to[i-1];
        to[i] = from[i] + tasks[i];
    }
}

