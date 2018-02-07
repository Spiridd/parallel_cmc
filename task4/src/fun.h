#ifndef FUN_H
#define FUN_H

void write_array_to_file(double *A, int m, int n, const char *filename);
void get_random_array(double *A, int N);
void* safe_malloc(size_t nbytes);

#endif

