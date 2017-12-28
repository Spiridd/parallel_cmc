// parallel implementation of the seive of Eratosthenes
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <mpi.h>

// maps index to number
int get_number(int i)
{
    assert(i>=0);
    return 3*i+5-(i%2);
}

// maps number to index
int get_index(int n)
{
    assert(n%6==1 || n%6==5);
    return n/3-1;
}

// returns size of list that can contain all the mappable
// numbers from 5 up to N (N is not included)
int get_size_of_candidates(int N)
{
    assert(N>5);
    const int rem = (N-1) % 6;
    int last_n;
    if (rem==1 || rem==5)       last_n = N-1;
    else if (rem==0 || rem==2)  last_n = N-2;
    else                        last_n = N-rem;
    return (int) (last_n-1)/3;
}

void* safe_malloc(int n)
{
    void* ret = malloc(n);
    if (ret == NULL)
    {
        fprintf(stderr, "Cannot allocate memory\n");
        exit(1);
    }
    return ret;
}

// sieves prepared array by itself
// array[0] is mapped to 5! (and is prime)
void sieve(char* is_prime, int size, int N)
{
    int i = 0;
    int p = get_number(i);
    int p2 = p*p;
    while (p2 < N)
    {
        if (is_prime[i])
        {
            // sieve mulitplies of p from p2 to N-1
            // TODO use better step
            for(int m=p2; m<N; m+=2*p)
            {
                if (m%6==1 || m%6==5)
                {
                    // m exists in the list
                    // and it is not prime
                    is_prime[get_index(m)] = 0;
                }
            }
        }
        ++i;
        p = get_number(i);
        p2 = p*p;
    }
    return;
}

int* get_small_primes(int N, int* primes_size)
{
    // allocate excess memory
    int* small_primes = safe_malloc(4*sizeof(int));
    small_primes[0] = 2;
    small_primes[1] = 3;
    small_primes[2] = 5;
    small_primes[3] = 7;
    if (N < 3)
    {
        *primes_size = 0;
    }
    else if (N < 4)
    {
        *primes_size = 1;
    } 
    else if (N < 6)
    {
        *primes_size = 2;
    } 
    else if (N < 8)
    {
        *primes_size = 3;
    } 
    else if (N < 12)
    {
        *primes_size = 4;
    }
    else
    {
        fprintf(stderr, "incorrect function parameter N = %d\n", N);
        exit(1);
    }
    return small_primes;
}

// serial straight algorithm
// returns pointer to array of primes up to N-1
// and puts its size into primes_size variable
int* get_primes(int N, int* primes_size)
{
    if (N < 12) return get_small_primes(N, primes_size);
    // 2 and 3 and their mulitplies are excluded
    // mapping functions are described above
    const int size = get_size_of_candidates(N);
    char* is_prime = (char*) safe_malloc(size*sizeof(char));

    // 0 - not prime, 1 - prime
    memset(is_prime, 1, size*sizeof(char));

    // marking not prime numbers
    sieve(is_prime, size, N);

    int count = 2; // 2, 3
    for(int i=0; i<size; ++i)
    {
        if (is_prime[i])    ++count;
    }

    *primes_size = count;
    int* primes = safe_malloc((*primes_size)*sizeof(int));
    primes[0] = 2;
    primes[1] = 3;
    primes[2] = 5;
    primes[3] = 7;
    int primes_idx = 4;
    for(int i=2; i<size; ++i)
    {
        if (is_prime[i])
        {
            primes[primes_idx] = get_number(i);
            ++primes_idx;
        }
    }
    free(is_prime);

    return primes;
}

int get_nearest_up_index(int a)
{
    int from = a+1;
    while (from % 6 != 1 && from % 6 != 5)  ++from;    
    return get_index(from);
}

int get_nearest_down_index(int b)
{
    int to = b-1;
    while (to % 6 != 1 && to % 6 != 5)  --to;    
    return get_index(to);
}

void sieve_with_primes(int* primes, int primes_size, char* is_prime, int min_index, int max_index)
{
    const int min_value = get_number(min_index);
    const int max_value = get_number(max_index);
    for(int i=0; i<primes_size; ++i)
    {
        const int p = primes[i];
        const int rem = min_value % p;
        int start;
        if (rem == 0)
        {
            // case: p=k*min_value, k is integer
            // as 2p, 3p, 4p are not in the list apriori
            // we start from 5p
            start = (min_value==p ? 5*p : min_value);
        }
        else
        {
            // case: p<min_value or p>min_value
            // as p must remain in the list, we start from 5*p
            start = (min_value<p ? 5*p : min_value-rem+p);
        }
        // when we start from 2p, 4p, and so on
        if (start % 2 == 0)
        {
            start += p;
        }
        for(int num=start; num<=max_value; num+=2*p)
        {
            if (num%6==1 || num%6==5)
            {
                // num is in the list
                // and is not prime
                is_prime[get_index(num)-min_index] = 0;
            }
        }
    }
    return;
}

// IMPORTANT it is service function
// it does not work when a<3 
// as it invokes some unnecessary difficulties
// returns primes in range (a; b) and sets size in primes_size
// and aux_primes through parameters
int* get_primes_in_between(int a, int b, int* primes_size,
       int** aux_primes, int* aux_primes_size)
{
    assert(a>=3); // 2 and 3 are excluded
    assert(b-a>1); // not empty interval
    // 2 and 3 and their mulitplies are excluded
    // mapping functions are described above
    const int min_index = get_nearest_up_index(a);
    const int max_index = get_nearest_down_index(b);
    const int size = max_index-min_index+1;
    char* is_prime = safe_malloc(size*sizeof(char));
    memset(is_prime, 1, size);

    // auxillary primes generation
    // get_primes returns primes up to N-1. So we add 1
    *aux_primes = get_primes((int)sqrt(b-1)+1, aux_primes_size);
    sieve_with_primes(*aux_primes, *aux_primes_size, is_prime, min_index, max_index);
    // count and create array of ints
    int count = 0;
    for(int i=0; i<size; ++i)
    {
        if (is_prime[i])    ++count;
    }
    *primes_size = count;
    int* primes = safe_malloc((*primes_size)*sizeof(int));
    int primes_idx = 0;
    for(int i=0; i<size; ++i)
    {
        if (is_prime[i])
        {
            primes[primes_idx] = get_number(min_index+i);
            ++primes_idx;
        }
    }
    free(is_prime);
        
    return primes;
}

// generate filename for output
char* get_res_filename(int rank)
{
    char* str = malloc(20);
    strcpy(str, "out/res_");
    sprintf(str+8, "%d.txt", rank);
    return str;
}

// generate filename for data
char* get_data_filename(int n_proc)
{
    char* str = malloc(20);
    strcpy(str, "out/data_");
    sprintf(str+9, "%d.txt", n_proc);
    return str;
}

// prints primes to stdout
void parallel_sieve(int N)
{
    // get rank and size of processors
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    
    // create data file
    //if (rank == 0)
        //fclose(fopen(get_data_filename(n_proc), "w"));

    // start counting time here
    const clock_t begin = clock();

    // get problem size and divide it
    const int sqrt_N1 = (int) sqrt(N-1);
    const int s = N-sqrt_N1-1; // size of (sqrt_N1; N)
    const int ds = s/n_proc; // piece of work for 1 processor
    // get filename for output
    char* res_filename = get_res_filename(rank);
    FILE* res_file = fopen(res_filename, "w");
    free(res_filename);
    //get filename for data (to keep time and rank)
    char* data_filename = get_data_filename(n_proc);
    FILE* data_file = fopen(data_filename, "a");
    free(data_filename);

    // calculate interval (a; b)
    if (rank > 0)
    {
        const int a = sqrt_N1 + (rank-1)*ds;
        const int b = sqrt_N1 + rank*ds + 1;
        int* aux_primes;
        int aux_primes_size;
        int primes_size;
        int* primes = get_primes_in_between(a, b, &primes_size, &aux_primes, &aux_primes_size);
        free(aux_primes);

        // stop counting time here
        const clock_t end = clock();
        const double time = (double)(end-begin)/CLOCKS_PER_SEC;
        fprintf(data_file, "%d,%d,%f\n", rank, N, time);

        for(int i=0; i<primes_size; ++i)
            fprintf(res_file, "%d ", primes[i]);
        free(primes);
    }
    else
    {
        // rank 0 (master)
        const int a = sqrt_N1 + (n_proc-1)*ds;
        const int b = N;
        int* aux_primes;
        int aux_primes_size;
        int primes_size;
        int* primes = get_primes_in_between(a, b, &primes_size, &aux_primes, &aux_primes_size);

        // stop counting time here
        const clock_t end = clock();
        const double time = (double)(end-begin)/CLOCKS_PER_SEC;
        fprintf(data_file, "%d,%d,%f\n", rank, N, time);

        for(int i=0; i<aux_primes_size; ++i)
            fprintf(res_file, "%d ", aux_primes[i]);
        free(aux_primes);

        // print answers 
        for(int i=0; i<primes_size; ++i)
            fprintf(res_file, "%d ", primes[i]);
        free(primes);
    }
    fclose(res_file);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    if (argc < 2)
    {
        fprintf(stderr, "provide size argument\n");
        return 1;
    }
    const int N = (int) strtol(argv[1], NULL, 10);
    parallel_sieve(N); 
    
    MPI_Finalize();
    return 0;
}
