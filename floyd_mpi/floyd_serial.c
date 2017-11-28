#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <sys/time.h>

double ts;

double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double) t.tv_sec + (double) t.tv_usec * 1E-6;
}

static inline void infinitize(int n, int *g)
{
    for (int i = 0; i < n * n; ++i)
        if (g[i] == 0)
            g[i] = 10000;
}

static inline void deinfinitize(int n, int *g)
{
    for (int i = 0; i < n * n; ++i)
        if (g[i] == 10000)
            g[i] = 0;
}

int *shortest_paths_s(int n, int *restrict g)
{
    int *restrict gnew = (int *) calloc(n * n, sizeof(int));
    memcpy(gnew, g, n * n * sizeof(int));
    infinitize(n, gnew);
    ts = wtime();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++) {
                if (gnew[j * n + k] > gnew[j * n + i] + gnew[i * n + k])
                    gnew[j * n + k] = gnew[j * n + i] + gnew[i * n + k];
            }
    }
    ts = wtime() - ts;
    for (int i = 0; i < n * n; i += n + 1)
        gnew[i] = 0;
    deinfinitize(n, gnew);
    return gnew;
}

int *gen_graph(int n)
{
    int *g = calloc(n * n, sizeof(int));
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i)
            if (j > i)
                g[j * n + i] = (rand() % 100);
            else
                g[j * n + i] = -1;
    }
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i)
            if (i > j)
                g[j * n + i] = g[i * n + j];
        g[j * n + j] = 0;
    }
    return g;
}

void write_matrix(const char *fname, int n, int *a)
{
    FILE *fp = fopen(fname, "w+");
    if (fp == NULL) {
        fprintf(stderr, "Could not open output file: %s\n", fname);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            fprintf(fp, "%d ", a[j * n + i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

int main(int argc, char **argv)
{
    int n = (argc > 1) ? atoi(argv[1]) : 100;
    const char *genname = (argc > 2) ? argv[2] : NULL;
    const char *sername = (argc > 3) ? argv[3] : NULL;
    if (n < 1) {
        fprintf(stderr, "Invalid size of graph");
        exit(EXIT_FAILURE);
    }

    int *g = gen_graph(n);
    if (genname)
        write_matrix(genname, n, g);
    int *gs = shortest_paths_s(n, g);

    printf("n:     %d\n", n);
    printf("Time ser (sec): %.6f\n", ts);

    if (sername)
        write_matrix(sername, n, gs);

    free(g);
    free(gs);
    return 0;
}
