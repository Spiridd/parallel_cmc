void write_array_to_file(double *A, int m, int n, const char *filename);
void get_random_array(double *A, int N);
void* safe_malloc(size_t nbytes);
void distribute_tasks(int n_tasks, int n_proc, int **tasks_ptr, int **from_ptr, int **to_ptr);

