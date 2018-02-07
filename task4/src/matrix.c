#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <mpi.h>
#include "fun.h"
#include "matrix.h"
#define N_ITER 5

void read_matrix_size_from_file(const char *filename, int *m, int *n)
{
    FILE *file = fopen(filename, "rb");
    if (file==NULL)
    {
        fprintf(stderr, "file %s cannot be opened\n", filename);
        exit(1);
    }
    int type;
    fread(&type, sizeof(int), 1, file);
    assert(type==1); // expect double only
    fread(m, sizeof(int), 1, file);
    assert(*m>0);
    fread(n, sizeof(int), 1, file);
    assert(*n>0);
    fclose(file);
}

void read_sizes_from_files(const char *file_A, const char *file_b, int *m, int *n)
{
    /* A: m*n, b: n*1 */
    read_matrix_size_from_file(file_A, m, n);
    int mb, nb;
    read_matrix_size_from_file(file_b, &mb, &nb);
    assert(*n==mb);
    assert(nb==1);
}

void read_matrices_portions_from_files(const char *file_A, const char *file_b, int from, int to, int is_col_blocked, int m, int n, double **A_ptr, double **b_ptr)
{
    // A: m*n, b: n*1 
    FILE *f1 = fopen(file_A, "rb");
    FILE *f2 = fopen(file_b, "rb");
    if (f1==NULL || f2==NULL)
    {
        fprintf(stderr, "one or both of these files cannot be opened: %s %s\n",
                file_A, file_b);
        exit(EXIT_FAILURE);
    }
    const int n_tasks = to-from;
    fseek(f1, 3*sizeof(int), SEEK_SET);
    fseek(f2, 3*sizeof(int), SEEK_SET);
    if (!is_col_blocked) // row blocked
    {
        *A_ptr = safe_malloc(n_tasks*n*sizeof(double));
        double * const A = *A_ptr;
        fseek(f1, from*n*sizeof(double), SEEK_CUR);
        fread(A, sizeof(double), n_tasks*n, f1);

        *b_ptr = safe_malloc(n*sizeof(double));
        double * const b = *b_ptr;
        fread(b, sizeof(double), n, f2);
    }
    else // coloumn blocked
    {
        *A_ptr = safe_malloc(m*n_tasks*sizeof(double));
        double * const A = *A_ptr;
        fseek(f1, from*sizeof(double), SEEK_CUR);
        for(int r=0; r<m; ++r)
        {
            fread(A+r*n_tasks, sizeof(double), n_tasks, f1);
            fseek(f1, (n-n_tasks)*sizeof(double), SEEK_CUR);
        }

        *b_ptr = safe_malloc(n_tasks*sizeof(double));
        double * const b = *b_ptr;
        fseek(f2, from*sizeof(double), SEEK_CUR);
        fread(b, sizeof(double), n_tasks, f2);
    }

    fclose(f1);
    fclose(f2);
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

void parallel_multiply(const char *file_A, const char *file_b, const char *file_c)
{
    int m, n;
    read_sizes_from_files(file_A, file_b, &m, &n);
    int n_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    int *tasks, *from, *to;

    if (m >= n)
    {
        // row blocked (why?)
        distribute_tasks(m, n_proc, &tasks, &from, &to);
        parallel_multiply_row_blocked(file_A, file_b, file_c, m, n, tasks, from, to);
    }
    else
    {
        // coloumn_blocked
        distribute_tasks(n, n_proc, &tasks, &from, &to);
        parallel_multiply_coloumn_blocked(file_A, file_b, file_c, m, n, tasks, from, to);
    }

    free(tasks);
    free(from);
    free(to);
}

// row blocked multiplication
// process 0 gathers the result and prints it into a file
void parallel_multiply_row_blocked(const char *file_A, const char *file_b, const char *file_c, int m, int n, const int *tasks, const int *from, const int *to)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        int n_proc;
        MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
        printf("n_proc: %d\n"
               "matrix: %dx%d\n"
               "row blocked\n", n_proc, m, n);
    }

    double *A, *b;
    read_matrices_portions_from_files(file_A, file_b, from[rank], to[rank], 0, m, n, &A, &b);
    double * const c = safe_malloc(tasks[rank]*sizeof(double));
    double *recvbuf = NULL;
    if (rank == 0)  recvbuf = safe_malloc(m*sizeof(double));
    clock_t reduce_clock = 0;
    //!TEST MPI_Wtime() vs clock()
    //const clock_t begin = clock();
    const double begin = MPI_Wtime(); 
    for(int iter=0; iter<N_ITER; ++iter)
    {
        /* calculating portion */
        for(int i=0; i<tasks[rank]; ++i)
        {
            const int i_n = i*n;
            double ci = 0.0;
            for(int k=0; k<n; ++k)
                ci += A[i_n+k] * b[k];
            c[i] = ci;
        }
        /* master process gathers data: just copies portions */
        const clock_t reduce_begin = clock();
        MPI_Gatherv(c, tasks[rank], MPI_DOUBLE, recvbuf, tasks, from, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        const clock_t reduce_end = clock();
        reduce_clock += reduce_end-reduce_begin;
    }
    //const clock_t end = clock();
    const double end = MPI_Wtime();
    free(A);
    free(b);
    free(c);
    //const double mytime = (double)(end-begin)/CLOCKS_PER_SEC/N_ITER;
    const double mytime = (end-begin)/N_ITER;
    double max_time;
    MPI_Reduce(&mytime, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    const double my_reduce = (double)reduce_clock/CLOCKS_PER_SEC/N_ITER;
    double max_reduce;
    MPI_Reduce(&my_reduce, &max_reduce, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("MPI time: %f\n"
               "Reduce time: %f\n", max_time, max_reduce);
        write_array_to_file(recvbuf, m, 1, file_c);
        free(recvbuf);
    }
}

// coloumn blocked multiplication
// process 0 gathers the result and prints it into a file
void parallel_multiply_coloumn_blocked(const char *file_A, const char *file_b, const char *file_c, int m, int n, const int *tasks, const int *from, const int *to)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        int n_proc;
        MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
        printf("n_proc: %d\n"
               "matrix: %dx%d\n"
               "coloumn blocked\n", n_proc, m, n);
    }

    double *A, *b;
    read_matrices_portions_from_files(file_A, file_b, from[rank], to[rank], 1, m, n, &A, &b);
    double * const c = safe_malloc(m*sizeof(double));
    double *recvbuf = NULL;
    if (rank == 0)  recvbuf = safe_malloc(m*sizeof(double));
    clock_t reduce_clock = 0;
    const clock_t begin = clock();
    for(int iter=0; iter<N_ITER; ++iter)
    {
        const int n_tasks = tasks[rank];
        /* calculating portion */
        for(int i=0; i<m; ++i)
        {
            const int i_offset = i*n_tasks;
            double ci = 0.0;
            for(int k=0; k<n_tasks; ++k)
                ci += A[i_offset+k] * b[k];
            c[i] = ci;
        }
        /* master process gathers data: sum of all the results */
        const clock_t reduce_begin = clock();
        MPI_Reduce(c, recvbuf, m, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        const clock_t reduce_end = clock();
        reduce_clock += reduce_end-reduce_begin;
    }
    const clock_t end = clock();
    free(A);
    free(b);
    free(c);
    const double mytime = (double)(end-begin)/CLOCKS_PER_SEC/N_ITER;
    double max_time;
    MPI_Reduce(&mytime, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    const double my_reduce = (double)reduce_clock/CLOCKS_PER_SEC/N_ITER;
    double max_reduce;
    MPI_Reduce(&my_reduce, &max_reduce, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("MPI time: %f\n"
               "Reduce time: %f\n", max_time, max_reduce);
        write_array_to_file(recvbuf, m, 1, file_c);
        free(recvbuf);
    }
}

