#ifndef MATRIX_H
#define MATRIX_H

void read_matrix_size_from_file(const char *filename, int *m, int *n);
void read_sizes_from_files(const char *file_A, const char *file_b, int *m, int *n);
void read_matrices_portions_from_files(const char *file_A, const char *file_b, int from, int to, int is_col_blocked, int m, int n, double **A_ptr, double **b_ptr);
void distribute_tasks(int n_tasks, int n_proc, int **tasks_ptr, int **from_ptr, int **to_ptr);
void parallel_multiply(const char *file_A, const char *file_b, const char *file_c);

// row blocked multiplication
// process 0 gathers the result and prints it into a file
void parallel_multiply_row_blocked(const char *file_A, const char *file_b, const char *file_c, int m, int n, const int *tasks, const int *from, const int *to);

// coloumn blocked multiplication
// process 0 gathers the result and prints it into a file
void parallel_multiply_coloumn_blocked(const char *file_A, const char *file_b, const char *file_c, int m, int n, const int *tasks, const int *from, const int *to);

#endif

