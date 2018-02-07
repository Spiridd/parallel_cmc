/*
 * generates mapping file for Blue Gene/P
 * for VN mode with 128 nodes (512 proc)
 * or for DUAL mode with 256 nodes (512 proc)
 * and writes it in my.map
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

void get_coords(int n, int *dims, int *coords)
{
    assert(n>=0);
    assert(n<512);
    const int cube_size = 512/dims[0];
    const int plane_size = cube_size/dims[1];
    const int line_size = plane_size/dims[2];
    assert(line_size==dims[3]);
    coords[0] = n/cube_size;
    coords[1] = (n%cube_size)/plane_size;
    coords[2] = ((n%cube_size)%plane_size)/line_size;
    coords[3] = n%line_size;
}

void make_map(int *dims, const char *filename)
{
    int numbers[512];
    for(int i=0; i<512; ++i)
        numbers[i] = i;
    // shuffle
    for(int i=0; i<512; ++i)
    {
        const int pos = rand() % 512;
        // swap i and pos
        const int temp = numbers[i];
        numbers[i] = numbers[pos];
        numbers[pos] = temp;
    }
    FILE *f = fopen(filename, "w");
    if (f == NULL)
    {
        fprintf(stderr, "can't open/create %s\n", filename);
        exit(EXIT_FAILURE);
    }
    // write to the file
    for(int i=0; i<512; ++i)
    {
        int coords[4];
        get_coords(numbers[i], dims, coords);
        fprintf(f, "%d %d %d %d\n", coords[0], coords[1], coords[2], coords[3]);
    }
    fclose(f);
}

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s mode filename\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    int dims[4];
    if (strcmp(argv[1], "VN") == 0)
    {
        dims[0] = 4;
        dims[1] = 4;
        dims[2] = 8;
        dims[3] = 4;
    }
    else if (strcmp(argv[1], "DUAL") == 0)
    {
        fprintf(stderr, "DUAL mode doesn't work for some reason\n");
        dims[0] = 8;
        dims[1] = 4;
        dims[2] = 8;
        dims[3] = 2;
    }
    else
    {
        fprintf(stderr, "incorrect mode: %s\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    
    make_map(dims, argv[2]);

    return 0;
}

