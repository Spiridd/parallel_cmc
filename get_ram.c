#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

int main()
{
    int f;
    if ((f = open("/proc/meminfo", O_RDONLY)) == -1)
    {
        fprintf(stderr, "error on opening\n");
        return 1;
    }

    char buf[BUFSIZ];
    int n;
    while ((n = read(f, buf, BUFSIZ)) > 0)
    {
        if (write(1, buf, n) != n)
        {
            fprintf(stderr, "error on output\n");
            return 1;
        }
    }
    close(f);

    return 0;
}
