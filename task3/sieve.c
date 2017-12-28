// serial implementation of the seive of Eratosthenes
// there 2 functions: naive and better one
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

// prints number of primes in range [2, N)
// returns their count
unsigned int get_count_of_primes(int N)
{
    // value = its index in is_prime array
    int * const is_prime = malloc(N * sizeof(int));
    if (is_prime == NULL)
    {
        fprintf(stderr, "Cannot allocate memory\n");
        exit(1);
    }

    // 0 - not prime, 1 - prime
    memset(is_prime, 1, N*sizeof(int));
    is_prime[0] = is_prime[1] = 0;

    // remove multiples of 2
    for(int i=4; i<N; i+=2)
        is_prime[i] = 0;
     
    // remove multiples of p
    int p2;
    for(int p=3; p2=p*p, p2<N; p+=2)
    {
        if (is_prime[p])
        {
            for(int i=p2; i<N; i+=2*p)
                is_prime[i] = 0;
        }
    }

    unsigned int count = 0;
    for(int p=0; p<N; ++p)
    {
        if (is_prime[p])
        {
            //printf("%d ", p);
            ++count;
        }
    }
    printf("\n");

    free(is_prime);

    return count;
}

unsigned int get_count_of_primes_naive(int N)
{
    int * is_prime = malloc(N * sizeof(int));
    if (is_prime == NULL)
    {
        fprintf(stderr, "Cannot allocate memory\n");
        exit(1);
    }

    for(int i=2; i<N; ++i)
        is_prime[i] = 1;
    is_prime[0] = is_prime[1] = 0;
    // remove multiples of 2
    for(int i=4; i<N; i+=2)
        is_prime[i] = 0;
    // remove multiples of p
    for(int p=3; p<N; p+=2)
    {
        if (is_prime[p])
        {
            for(int i=2*p; i<N; i+=p)
                is_prime[i] = 0;
        }
    }

    unsigned int count = 0;
    for(int p=0; p<N; ++p)
        if (is_prime[p])        ++count;

    free(is_prime);

    return count;
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "provide size argument\n");
        return 1;
    }
    const int N = (int) strtol(argv[1], NULL, 10);
    if (N < 2)
    {
        fprintf(stderr, "size must be positive integer bigger than 1\n");
        return 1;
    }
    // speed comparison
    static const int n_iter = 1;

    clock_t begin, end;
    int count;
    
    begin = clock();
    for(int i=0; i<n_iter; ++i)
        count = get_count_of_primes_naive(N);
    end = clock();
    const double time_naive = (double)(end-begin)/CLOCKS_PER_SEC/n_iter;
    printf("naive time = %f\n", time_naive);
    printf("naive result = %d\n", count);

    begin = clock();
    for(int i=0; i<n_iter; ++i)
        count = get_count_of_primes(N);
    end = clock();
    const double time = (double)(end-begin)/CLOCKS_PER_SEC/n_iter;
    printf("smart time = %f\n", time);
    printf("smart result = %d\n", count);

    printf("----------\n"
            "speedup: %f\n", time_naive/time);
    return 0;
}
