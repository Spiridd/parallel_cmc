#ifndef MATRIX_H
#define MATRIX_H

void distribute_tasks(int n_tasks, int n_proc, int **tasks_ptr, int **from_ptr, int **to_ptr);
void* safe_malloc(size_t nbytes);
void get_cartesian_comm(MPI_Comm *comm);
int get_matrix_size_from_file(MPI_Comm comm, const char *filename);
void read_matrices_from_files(MPI_Comm cart_comm, const char *filename_A, const char *filename_B, double **A_ptr, double **B_ptr, int *i0, int *i1, int *j0, int *j1, int *k0, int *k1, int n);
void parallel_multiply(const char *filename_A, const char *filename_B, const char *filename_C);

#endif

