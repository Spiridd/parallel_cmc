// each matrix is square;
// only equal-sized matrix can be multiplied
#include <stdlib.h>
#include "matrix.h"

/*******************
  general functions
********************/

void get_random_array(float* A, int n)
{
    for(int i=0; i<n; ++i)
    {
        A[i] = (float) rand()/RAND_MAX;
    }
}

void set_zeros(float* A, int n)
{
    for(int i=0; i<n; ++i)
    {
        A[i] = 0.0;
    }
}

/*******************
  naive algorithms
********************/

void dgemm_ijk(float* A, float* B, float* C, int n)
{
    for(int i=0; i<n; ++i)
    {
        const int i_n = i*n;
        for(int j=0; j<n; ++j)
        {
            float cij = 0.0;
            for(int k=0; k<n; ++k)
                cij += A[i_n+k] * B[k*n+j];
            C[i_n+j] = cij;
        }
    }
}

void dgemm_jik(float* A, float* B, float* C, int n)
{
    for(int j=0; j<n; ++j)
    {
        for(int i=0; i<n; ++i)
        {
            const int i_n = i*n;
            float cij = 0.0;
            for(int k=0; k<n; ++k)
                cij += A[i_n+k] * B[k*n+j];
            C[i_n+j] = cij;
        }
    }
}

void dgemm_ikj(float* A, float* B, float* C, int n)
{
    for(int i=0; i<n; ++i)
    {
        const int i_n = i*n;
        for(int k=0; k<n; ++k)
        {
            const int k_n = k*n;
            const float a = A[i_n+k];
            for(int j=0; j<n; ++j)
                C[i_n+j] += a * B[k_n+j];
        }
    }
}

void dgemm_kij(float* A, float* B, float* C, int n)
{
    for(int k=0; k<n; ++k)
    {
        const int k_n = k*n;
        for(int i=0; i<n; ++i)
        {
            const int i_n = i*n;
            const float a = A[i_n+k];
            for(int j=0; j<n; ++j)
            {
                C[i_n+j] += a * B[k_n+j];
            }
        }
    }
}

void dgemm_kji(float* A, float* B, float* C, int n)
{
    for(int k=0; k<n; ++k)
    {
        const int k_n = k*n;
        for(int j=0; j<n; ++j)
        {
            const float b = B[k_n+j];
            for(int i=0; i<n; ++i)
                C[i*n+j] += A[i*n+k] * b; 
        }
    }
}

void dgemm_jki(float* A, float* B, float* C, int n)
{
    for(int j=0; j<n; ++j)
    {
        for(int k=0; k<n; ++k)
        {
            const int k_n = k*n;
            const float b = B[k_n+j];
            for(int i=0; i<n; ++i)
                C[i*n+j] += A[i*n+k] * b;
        }
    }
}

/*******************
  square algorithms (L1 blocked)
********************/

void square_dgemm_ijk_ijk(float* A, float* B, float* C, int n, int bsz)
{
    for(int i=0; i<n; i+=bsz)
    {
        const int M = min(bsz, n-i);
        for(int j=0; j<n; j+=bsz)
        {
            const int N = min(bsz, n-j);
            for(int k=0; k<n; k+=bsz)
            {
                // correcting case of the last block
                const int K = min(bsz, n-k);
                // block multiplication ikj (fastest)
                for(int i1=0; i1<M; ++i1)
                {
                    const int i_offset = (i+i1)*n;
                    for(int k1=0; k1<K; ++k1)
                    {
                        const float a = A[i_offset+k+k1];
                        for(int j1=0; j1<N; ++j1)
                            C[i_offset+j+j1] += a * B[(k+k1)*n+j+j1];
                    }
                }
            }
        }
    }
}

void square_dgemm_ijk_ijk_opt(float* A, float* B, float* C, int n, int bsz)
{
    for(int i=0; i<n; i+=bsz)
    {
        for(int j=0; j<n; j+=bsz)
            for(int k=0; k<n; k+=bsz)
            {
                // correcting case of the last block
                const int M = min(bsz, n-i);
                const int N = min(bsz, n-j);
                const int K = min(bsz, n-k);
                // block multiplication
                for(int i1=0; i1<M; ++i1)
                {
                    const int i_offset = (i+i1)*n;
                    for(int j1=0; j1<N; ++j1)
                    {
                        float cij = 0.0;
                        int k1 = 0;
                        for(; k1<K-8; k1+=8)
                        {
                            const float d0 = A[i_offset+k+k1+0] * B[(k+k1+0)*n+j+j1];
                            const float d1 = A[i_offset+k+k1+1] * B[(k+k1+1)*n+j+j1];
                            const float d2 = A[i_offset+k+k1+2] * B[(k+k1+2)*n+j+j1];
                            const float d3 = A[i_offset+k+k1+3] * B[(k+k1+3)*n+j+j1];
                            const float d4 = A[i_offset+k+k1+4] * B[(k+k1+4)*n+j+j1];
                            const float d5 = A[i_offset+k+k1+5] * B[(k+k1+5)*n+j+j1];
                            const float d6 = A[i_offset+k+k1+6] * B[(k+k1+6)*n+j+j1];
                            const float d7 = A[i_offset+k+k1+7] * B[(k+k1+7)*n+j+j1];
                            cij += (d0+d1+d2+d3+d4+d5+d6+d7);
                        }

                        for( ; k1<K; ++k1)
                            cij += A[i_offset+k+k1] * B[(k+k1)*n+j+j1];

                        C[i_offset+j+j1] += cij;
                    }
                }
            }
    }
}

void square_dgemm_ijk_with_call_ijk_opt(float* A, float* B, float* C,
                                        int n, int bsz)
{
    for(int i=0; i<n; i+=bsz)
    {
        const int i_offset = i*n;
        const int M = min(bsz, n-i);
        for(int j=0; j<n; j+=bsz)
        {
            const int N = min(bsz, n-j);
            for(int k=0; k<n; k+=bsz)
            {
                const int K = min(bsz, n-k);
                calculate_block_ijk_opt(A+i_offset+k, B+k*n+j, C+i_offset+j,
                                        M, N, K, n);
            }
        }
    }
}

void square_dgemm_ikj_with_call_ijk_opt(float* A, float* B, float* C,
                                        int n, int bsz)
{
    for(int i=0; i<n; i+=bsz)
    {
        const int i_offset = i*n;
        const int M = min(bsz, n-i);
        for(int k=0; k<n; k+=bsz)
        {
            const int K = min(bsz, n-k);
            for(int j=0; j<n; j+=bsz)
            {
                const int N = min(bsz, n-j);
                calculate_block_ijk_opt(A+i_offset+k, B+k*n+j, C+i_offset+j,
                                        M, N, K, n);
            }
        }
    }
}

void square_dgemm_ikj_with_call_ikj(float* A, float* B, float* C,
                                        int n, int bsz)
{
    for(int i=0; i<n; i+=bsz)
    {
        const int i_offset = i*n;
        const int M = min(bsz, n-i);
        for(int k=0; k<n; k+=bsz)
        {
            const int K = min(bsz, n-k);
            for(int j=0; j<n; j+=bsz)
            {
                const int N = min(bsz, n-j);
                calculate_block_ikj(A+i_offset+k, B+k*n+j, C+i_offset+j,
                                        M, N, K, n);
            }
        }
    }
}

void square_dgemm_ijk_with_call_ikj(float* A, float* B, float* C,
                                        int n, int bsz)
{
    for(int i=0; i<n; i+=bsz)
    {
        const int i_offset = i*n;
        const int M = min(bsz, n-i);
        for(int j=0; j<n; j+=bsz)
        {
            const int N = min(bsz, n-j);
            for(int k=0; k<n; k+=bsz)
            {
                const int K = min(bsz, n-k);
                calculate_block_ikj(A+i_offset+k, B+k*n+j, C+i_offset+j,
                                        M, N, K, n);
            }
        }
    }
}

void square_dgemm_ijk_with_call_ijk(float* A, float* B, float* C,
                                        int n, int bsz)
{
    for(int i=0; i<n; i+=bsz)
    {
        const int i_offset = i*n;
        const int M = min(bsz, n-i);
        for(int j=0; j<n; j+=bsz)
        {
            const int N = min(bsz, n-j);
            for(int k=0; k<n; k+=bsz)
            {
                const int K = min(bsz, n-k);
                calculate_block_ijk(A+i_offset+k, B+k*n+j, C+i_offset+j,
                                        M, N, K, n);
            }
        }
    }
}

/*******************
  square algorithms: block functions (L1 blocked)
********************/

static void calculate_block_ijk(float* A, float* B, float* C,
                                    int M, int N, int K, int n)
{
    for(int i=0; i<M; ++i)
    {
        const int i_offset = i*n;
        for(int j=0; j<N; ++j)
        {
            float cij = 0.0;
            for(int k=0; k<K; ++k)
                cij += A[i_offset+k] * B[k*n+j];
            C[i_offset+j] += cij;
        }
    }
}

static void calculate_block_ikj(float* A, float* B, float* C,
                                    int M, int N, int K, int n)
{
    for(int i=0; i<M; ++i)
    {
        const int i_offset = i*n;
        for(int k=0; k<K; ++k)
        {
            const int k_offset = k*n;
            const double a = A[i_offset+k];
            for(int j=0; j<N; ++j)
                C[i_offset+j] += a * B[k_offset+j];
        }
    }
}

static void calculate_block_ijk_opt(float* A, float* B, float* C,
                                    int M, int N, int K, int n)
{
    for(int i=0; i<M; ++i)
    {
        const int i_offset = i*n;
        for(int j=0; j<N; ++j)
        {
            float cij = 0.0;
            int k = 0; 
            for( ; k<K-8; k+=8)
            {
                const float d0 = A[i_offset+k+0] * B[(k+0)*n+j];
                const float d1 = A[i_offset+k+1] * B[(k+1)*n+j];
                const float d2 = A[i_offset+k+2] * B[(k+2)*n+j];
                const float d3 = A[i_offset+k+3] * B[(k+3)*n+j];
                const float d4 = A[i_offset+k+4] * B[(k+4)*n+j];
                const float d5 = A[i_offset+k+5] * B[(k+5)*n+j];
                const float d6 = A[i_offset+k+6] * B[(k+6)*n+j];
                const float d7 = A[i_offset+k+7] * B[(k+7)*n+j];
                cij += (d0+d1+d2+d3+d4+d5+d6+d7);
            }
            
            for( ; k<K; ++k)
                cij += A[i_offset+k] * B[k*n+j];

            C[i_offset+j] += cij;
        }
    }
}

