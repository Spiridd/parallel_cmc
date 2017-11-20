// the seive of Eratosthenes
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

unsigned int get_count_of_primes(long long N)
{
    int* const is_prime = malloc(N * sizeof(int));

    for(long long i=2; i<N; ++i)
        is_prime[i] = 1;
    is_prime[0] = is_prime[1] = 0;
    // remove powers of 2
    for(long long i=4; i<N; i+=2)
        is_prime[i] = 0;
    
    // find sqrt(N)
    int sqrt_N = 1;
    while (sqrt_N*sqrt_N < N)   ++sqrt_N;
    // remove divided by p
    for(long long p=3; p<=sqrt_N; p+=2)
    {
        for(long long i=p*p; i<N; i+=2*p)
            is_prime[i] = 0;
    }

    unsigned int count = 0;
    for(long long p=0; p<N; ++p)
        if (is_prime[p])        ++count;

    free(is_prime);

    return count;
}

unsigned int get_count_of_primes_naive(long long N)
{
    int* is_prime = malloc(N * sizeof(int));

    for(long long i=2; i<N; ++i)
        is_prime[i] = 1;
    is_prime[0] = is_prime[1] = 0;
    // remove powers of 2
    for(long long i=4; i<N; i+=2)
        is_prime[i] = 0;
    // remove divided by p
    for(long long p=3; p<N; p+=2)
    {
        for(long long i=2*p; i<N; i+=p)
            is_prime[i] = 0;
    }

    unsigned int count = 0;
    for(long long p=0; p<N; ++p)
        if (is_prime[p])        ++count;

    free(is_prime);

    return count;
}

int main(int argc, char** argv)
{
    const long long N = (long long) strtol(argv[1], NULL, 10);
    static const int n_iter = 1;

    clock_t begin, end;
    int count;
    
    begin = clock();
    for(int i=0; i<n_iter; ++i)
        count = get_count_of_primes_naive(N);
    end = clock();
    const double time_naive = (double)(end-begin)/CLOCKS_PER_SEC/n_iter;
    printf("naive time = %lf\n", time_naive);
    printf("naive result = %d\n", count);

    begin = clock();
    for(int i=0; i<n_iter; ++i)
        count = get_count_of_primes(N);
    end = clock();
    const double time = (double)(end-begin)/CLOCKS_PER_SEC/n_iter;
    printf("smart time = %lf\n", time);
    printf("smart result = %d\n", count);

    return 0;
}
