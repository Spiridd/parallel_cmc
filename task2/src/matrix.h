#ifndef H_MATRIX
#define H_MATRIX

#define min(a, b) (a<=b?a:b)

/*******************
  general functions
********************/

void get_random_array(float* A, int n);

void set_zeros(float* A, int n);

/*******************
  naive algorithms
********************/

void dgemm_ijk(float* A, float* B, float* C, int n);
void dgemm_jik(float* A, float* B, float* C, int n);
void dgemm_ikj(float* A, float* B, float* C, int n);
void dgemm_kij(float* A, float* B, float* C, int n);
void dgemm_kji(float* A, float* B, float* C, int n);
void dgemm_jki(float* A, float* B, float* C, int n);

/*******************
  square algorithms (L1 blocked)
********************/

void square_dgemm_ijk_ijk(float* A, float* B, float* C, int n, int bsz);
void square_dgemm_ijk_ijk_opt(float* A, float* B, float* C, int n, int bsz);

void square_dgemm_ijk_with_call_ijk(float* A, float* B, float* C,
                                        int n, int bsz);
void square_dgemm_ijk_with_call_ijk_opt(float* A, float* B, float* C,
                                        int n, int bsz);
void square_dgemm_ijk_with_call_ikj(float* A, float* B, float* C,
                                        int n, int bsz);
void square_dgemm_ikj_with_call_ikj(float* A, float* B, float* C,
                                        int n, int bsz);
void square_dgemm_ikj_with_call_ijk_opt(float* A, float* B, float* C,
                                        int n, int bsz);
/*******************
  square algorithms: block functions (L1 blocked)
********************/

static void calculate_block_ijk(float* A, float* B, float* C,
                                    int M, int N, int K, int n);
static void calculate_block_ikj(float* A, float* B, float* C,
                                    int M, int N, int K, int n);
static void calculate_block_ijk_opt(float* A, float* B, float* C,
                                    int M, int N, int K, int n);

#endif

