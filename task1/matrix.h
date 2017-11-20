void write_time_to_binary_2D(double times[], size_t count, char *filename);

void write_matrix_to_binary_2D(double **A, int n, char *filename);

void print_matrix_2D(double **A, int n);

void get_random_matrix_2D(double **M, int n);

void set_zeros_2D(double **A, int n);

double** allocate_matrix_2D(int n);

void free_matrix_2D(double **A, int n);

void multiply_ijk_2D(double **A, double **B, double **C, int n);
void multiply_jik_2D(double **A, double **B, double **C, int n);
void multiply_ikj_2D(double **A, double **B, double **C, int n);
void multiply_kij_2D(double **A, double **B, double **C, int n);
void multiply_kji_2D(double **A, double **B, double **C, int n);
void multiply_jki_2D(double **A, double **B, double **C, int n);

double measure_time_2D(void (*multiplier)(double**, double**, double**, int), double **A, double **B, double **C, int n, int n_iter);
