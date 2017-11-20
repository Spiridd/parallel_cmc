#ifndef H_MATRIX
#define H_MATRIX
#define double float

void write_time_to_binary_2D(double*, size_t, char*);

void write_matrix_to_binary_2D(double**, int, char*);

void print_matrix_2D(double**, int);

void get_random_matrix_2D(double**, int);
void get_random_matrix_1D(double*, int);

void set_zeros_2D(double**, int);
void set_zeros_1D(double*, int);

double** allocate_matrix_2D(int);
double* allocate_matrix_1D(int);

void free_matrix_2D(double**, int);
void free_matrix_1D(double*, int);

void multiply_ijk_2D(double**, double**, double**, int);
void multiply_jik_2D(double**, double**, double**, int);
void multiply_ikj_2D(double**, double**, double**, int);
void multiply_kij_2D(double**, double**, double**, int);
void multiply_kji_2D(double**, double**, double**, int);
void multiply_jki_2D(double**, double**, double**, int);

void multiply_ijk_1D(double*, double*, double*, int);
void multiply_jik_1D(double*, double*, double*, int);
void multiply_ikj_1D(double*, double*, double*, int);
void multiply_kij_1D(double*, double*, double*, int);
void multiply_kji_1D(double*, double*, double*, int);
void multiply_jki_1D(double*, double*, double*, int);

double measure_time_2D(void (*multiplier)(double**, double**, double**, int), double**, double**, double**, int, int);
double measure_time_1D(void (*multiplier)(double*, double*, double*, int), double*, double*, double*, int, int);

void calculate_block_1D(double*, double*, double*, int, int, int, int);
void multiply_block_1D(double*, double*, double*, int, int);
void square_dgemm_1D(double*, double*, double*, int, int);

#endif
