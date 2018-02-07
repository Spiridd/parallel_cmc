#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <mpi.h>
#include "matrix.h"

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

void get_cartesian_comm(MPI_Comm *comm)
{
    int n_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int dims[] = {0, 0, 0};
    //! this function doesn't work properly
    //MPI_Dims_create(n_proc, 3, dims);
    //!CRUTCH
    dims[0] = dims[1] = dims[2] = 2;
    if (dims[0]*dims[1]*dims[2] != n_proc)
    {
        if (rank==0)    fprintf(stderr, "cannot use %d processors, only %d (crutch)\n", n_proc, dims[0]*dims[1]*dims[2]);
        exit(EXIT_FAILURE);
    }
    if (dims[0]!=dims[1] || dims[0]!=dims[2])
    {
        if (rank == 0)    fprintf(stderr, "Number of processes must be cube\n");
        exit(EXIT_FAILURE);
    }
    int periods[] = {0, 0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, comm);
    if (*comm == MPI_COMM_NULL)
    {
        if (rank == 0)  fprintf(stderr, "cannot create 3D cartesian topology with %d processes\n", n_proc);
        exit(EXIT_FAILURE);
    }
    if (rank==0)    printf("Using %d processors\n", n_proc);
}

int get_matrix_size_from_file(MPI_Comm comm, const char *filename)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    MPI_File fh;
    if (MPI_File_open(comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
    {
        if (rank==0)    fprintf(stderr, "Cannot open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
    int data[3];
    MPI_File_read(fh, data, 3, MPI_INT, MPI_STATUS_IGNORE);
    if (data[0] != 1)
    {
        if (rank==0)    fprintf(stderr, "matrix type must be double\n");
        exit(EXIT_FAILURE);
    }
    if (data[1] != data[2])
    {
        if (rank==0)    fprintf(stderr, "matrix must be square\n");
        exit(EXIT_FAILURE);
    }
    MPI_File_close(&fh);
    return data[1];
}

void read_all_matrices_from_files(MPI_Comm cart_comm, const char *filename_A, const char *filename_B, double **A_ptr, double **B_ptr, int *i0_ptr, int *i1_ptr, int *j0_ptr, int *j1_ptr, int *k0_ptr, int *k1_ptr, int n)
{
    // collective read with derived datatypes
    // get distribution
    int rank;
    MPI_Comm_rank(cart_comm, &rank);
    int coords[3]; // IJK coordinates of process
    int dims[3], periods[3];
    MPI_Cart_get(cart_comm, 3, dims, periods, coords);
    int *tasks, *from, *to;
    distribute_tasks(n, dims[0], &tasks, &from, &to);
    free(tasks);

    const int i0 = *i0_ptr = from[coords[0]];
    const int j0 = *j0_ptr = from[coords[1]];
    const int k0 = *k0_ptr = from[coords[2]];
    const int i1 = *i1_ptr = to[coords[0]];
    const int j1 = *j1_ptr = to[coords[1]];
    const int k1 = *k1_ptr = to[coords[2]];

    free(from);
    free(to);

    *A_ptr = safe_malloc((i1-i0)*(k1-k0)*sizeof(double));
    double * const A = *A_ptr;
    *B_ptr = safe_malloc((j1-j0)*(k1-k0)*sizeof(double));
    double * const B = *B_ptr;

    // read the A portion
    MPI_File fh_A;
    if (MPI_File_open(cart_comm, filename_A, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_A))
    {
        if (rank==0)    fprintf(stderr, "Cannot open file %s\n", filename_A);
        exit(EXIT_FAILURE);
    }
    const int width_A = k1-k0;
    MPI_Datatype matrix_A_type;
    const int row_A_count = i1-i0;
    int * const displs_A = safe_malloc(row_A_count*sizeof(int));
    displs_A[0] = i0*n+k0;
    for(int i=1; i<row_A_count; ++i)  displs_A[i] = displs_A[i-1]+n;
    if (MPI_Type_create_indexed_block(row_A_count, width_A, displs_A, MPI_DOUBLE, &matrix_A_type))
    {
        if (rank==0)    fprintf(stderr, "cannot create type\n");
        exit(EXIT_FAILURE);
    }
    free(displs_A);
    MPI_Type_commit(&matrix_A_type);
    MPI_File_set_view(fh_A, 3*sizeof(int), MPI_DOUBLE, matrix_A_type, "native", MPI_INFO_NULL);
    MPI_File_read_at_all(fh_A, 0, A, row_A_count*width_A, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_Type_free(&matrix_A_type);
    MPI_File_close(&fh_A);
    
    // read the B portion
    MPI_File fh_B;
    if (MPI_File_open(cart_comm, filename_B, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_B))
    {
        if (rank==0)    fprintf(stderr, "Cannot open file %s\n", filename_B);
        exit(EXIT_FAILURE);
    }
    const int width_B = j1-j0;
    MPI_Datatype matrix_B_type;
    const int row_B_count = width_A;
    int * const displs_B = safe_malloc(row_B_count*sizeof(int));
    displs_B[0] = k0*n+j0;
    for(int k=1; k<row_B_count; ++k)  displs_B[k] = displs_B[k-1]+n;
    if (MPI_Type_create_indexed_block(row_B_count, width_B, displs_B, MPI_DOUBLE, &matrix_B_type))
    {
        if (rank==0)    fprintf(stderr, "cannot create type\n");
        exit(EXIT_FAILURE);
    }
    free(displs_B);
    MPI_Type_commit(&matrix_B_type);
    MPI_File_set_view(fh_B, 3*sizeof(int), MPI_DOUBLE, matrix_B_type, "native", MPI_INFO_NULL);
    MPI_File_read_at_all(fh_B, 0, B, row_B_count*width_B, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_Type_free(&matrix_B_type);
    MPI_File_close(&fh_B);
}

void read_matrices_from_files(MPI_Comm cart_comm, const char *filename_A, const char *filename_B, double **A_ptr, double **B_ptr, int *i0_ptr, int *i1_ptr, int *j0_ptr, int *j1_ptr, int *k0_ptr, int *k1_ptr, int n)
{
    // reads matrices with non-collective reads, row by row
    // get distribution
    int rank;
    MPI_Comm_rank(cart_comm, &rank);
    int coords[3]; // IJK coordinates of process
    int dims[3], periods[3];
    MPI_Cart_get(cart_comm, 3, dims, periods, coords);
    int *tasks, *from, *to;
    distribute_tasks(n, dims[0], &tasks, &from, &to);
    free(tasks);

    const int i0 = *i0_ptr = from[coords[0]];
    const int j0 = *j0_ptr = from[coords[1]];
    const int k0 = *k0_ptr = from[coords[2]];
    const int i1 = *i1_ptr = to[coords[0]];
    const int j1 = *j1_ptr = to[coords[1]];
    const int k1 = *k1_ptr = to[coords[2]];

    free(from);
    free(to);

    *A_ptr = safe_malloc((i1-i0)*(k1-k0)*sizeof(double));
    double * const A = *A_ptr;
    *B_ptr = safe_malloc((j1-j0)*(k1-k0)*sizeof(double));
    double * const B = *B_ptr;

    // read the A portion
    MPI_File fh_A;
    if (MPI_File_open(cart_comm, filename_A, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_A))
    {
        if (rank==0)    fprintf(stderr, "Cannot open file %s\n", filename_A);
        exit(EXIT_FAILURE);
    }
    const int width_A = k1-k0;
    MPI_File_set_view(fh_A, 3*sizeof(int), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
    for(int i=i0; i<i1; ++i)
    {
        MPI_File_read_at(fh_A, i*n+k0, A+width_A*(i-i0), width_A, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }
    MPI_File_close(&fh_A);
    
    // read the B portion
    MPI_File fh_B;
    if (MPI_File_open(cart_comm, filename_B, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_B))
    {
        if (rank==0)    fprintf(stderr, "Cannot open file %s\n", filename_B);
        exit(EXIT_FAILURE);
    }
    const int width_B = j1-j0;
    MPI_File_set_view(fh_B, 3*sizeof(int), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
    for(int k=k0; k<k1; ++k)
    {
        MPI_File_read_at(fh_B, k*n+j0, B+width_B*(k-k0), width_B, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }
    MPI_File_close(&fh_B);
}

void write_matrix_in_file(const char *filename, double *C, int i0, int i1, int j0, int j1, int n, MPI_Comm comm)
{
    // each processor writes its portion row by row independently
    MPI_File fh;
    MPI_File_open(comm, filename, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    const static int d_type = 1; // see in fun.c
    int data[] = {d_type, n, n};
    MPI_File_write_at(fh, 0, data, 3, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_set_view(fh, 3*sizeof(int), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
    for(int i=i0; i<i1; ++i)
    {
        MPI_File_write_at(fh, i*n+j0, C+(i-i0)*(j1-j0), j1-j0, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }
    MPI_File_close(&fh);
}

void write_all_matrix_in_file(const char *filename, const double *C, int i0, int i1, int j0, int j1, int n, MPI_Comm comm)
{
    // processors perform collective write with derived datatypes
    MPI_File fh;
    MPI_File_open(comm, filename, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    const static int d_type = 1; // see in fun.c
    int data[] = {d_type, n, n};
    MPI_File_write_at_all(fh, 0, data, 3, MPI_INT, MPI_STATUS_IGNORE);
    const int width = j1-j0;
    MPI_Datatype matrix_type;
    const int row_count = i1-i0;
    int * const displs = safe_malloc(row_count*sizeof(int));
    displs[0] = i0*n+j0;
    for(int i=1; i<row_count; ++i)  displs[i] = displs[i-1]+n;
    if (MPI_Type_create_indexed_block(row_count, width, displs, MPI_DOUBLE, &matrix_type))
    {
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank==0)    fprintf(stderr, "cannot create type\n");
        exit(EXIT_FAILURE);
    }
    free(displs);
    MPI_Type_commit(&matrix_type);
    MPI_File_set_view(fh, 3*sizeof(int), MPI_DOUBLE, matrix_type, "native", MPI_INFO_NULL);
    MPI_File_write_at_all(fh, 0, C, row_count*width, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_Type_free(&matrix_type);
    MPI_File_close(&fh);
}

void dgemm_ikj(double *A, double *B, double *C, int n, int m, int h)
{
    for(int i=0; i<n; ++i)
    {
        const int i_m = i*m;
        const int i_h = i*h;
        for(int k=0; k<m; ++k)
        {
            const int k_h = k*h;
            const double a = A[i_m+k];
            for(int j=0; j<h; ++j)
                C[i_h+j] += a * B[k_h+j];
        }
    }
}

void parallel_multiply(const char *filename_A, const char *filename_B, const char *filename_C)
{
    /* Cannon's algorithm */
    // Create cartesian topology
    MPI_Comm cart_comm;
    get_cartesian_comm(&cart_comm);
    int cart_rank;
    MPI_Comm_rank(cart_comm, &cart_rank);
    // Read matrices in parallel
    const int n = get_matrix_size_from_file(cart_comm, filename_A);
    assert(n == get_matrix_size_from_file(cart_comm, filename_B));
    double *A, *B;
    int i0, i1, j0, j1, k0, k1;
    const clock_t begin1 = clock();
    read_all_matrices_from_files(cart_comm, filename_A, filename_B, &A, &B, &i0, &i1, &j0, &j1, &k0, &k1, n);
    const clock_t end1 = clock();
    const double time1 = (double)(end1-begin1)/CLOCKS_PER_SEC;
    double max_time;
    MPI_Reduce(&time1, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, cart_comm);
    if (cart_rank==0)   printf("MPI read time: %f\n", max_time);
    // Calculate in parallel (MPI_Cart_sub)
    const int size_C = (i1-i0)*(j1-j0);
    double * const C = safe_malloc(size_C*sizeof(double));
    for(int i=0; i<size_C; ++i) C[i] = 0.0;
    const clock_t begin = clock();
    dgemm_ikj(A, B, C, i1-i0, k1-k0, j1-j0);
    free(A);
    free(B);
    // reduce in third dimension
    int dim_subs[] = {0, 0, 1}; 
    MPI_Comm linear_comm;
    MPI_Cart_sub(cart_comm, dim_subs, &linear_comm);
    int rank;
    MPI_Comm_rank(linear_comm, &rank);
    double *global_C = NULL; 
    if (rank == 0) global_C = safe_malloc(size_C*sizeof(double)); 
    MPI_Reduce(C, global_C, size_C, MPI_DOUBLE, MPI_SUM, 0, linear_comm);
    const clock_t end = clock();
    const double time = (double)(end-begin)/CLOCKS_PER_SEC;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, cart_comm);
    if (cart_rank==0)   printf("MPI time: %f\n", max_time);
    free(C);
    MPI_Comm_free(&linear_comm);
    // Write matrices within a new communicator
    dim_subs[0] = dim_subs[1] = 1;
    dim_subs[2] = 0;
    MPI_Comm square_comm;
    MPI_Cart_sub(cart_comm, dim_subs, &square_comm);
    const clock_t begin2 = clock();
    if (rank == 0)  write_all_matrix_in_file(filename_C, global_C, i0, i1, j0, j1, n, square_comm);
    const clock_t end2 = clock();
    const double time2 = (double)(end2-begin2)/CLOCKS_PER_SEC;
    MPI_Reduce(&time2, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, cart_comm);
    if (cart_rank==0)   printf("MPI write time: %f\n", max_time);
    free(global_C);
    MPI_Comm_free(&square_comm);
    MPI_Comm_free(&cart_comm);
}

