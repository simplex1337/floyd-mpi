#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include <sys/time.h>

double tp;

double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double) t.tv_sec + (double) t.tv_usec * 1E-6;
}

void floyd(int *g, int n, int threads)
{
    for (int k = 0; k < n; ++k) {
#pragma omp parallel for num_threads(threads)
        for (int i = 0; i < n; ++i) {
            int v = g[i * n + k];
            for (int j = 0; j < n; ++j) {
                int val = v + g[k * n + j];
                if (g[i * n + j] > val) {
                    g[i * n + j] = val;
                }
            }
        }
    }
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

int *shortest_paths_p(int n, int *restrict g)
{
    int *restrict gnew = (int *) calloc(n * n, sizeof(int));
    memcpy(gnew, g, n * n * sizeof(int));
    infinitize(n, gnew);
    tp = omp_get_wtime();
    floyd(gnew, n, omp_get_max_threads());
    tp = omp_get_wtime() - tp;
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
    MPI_Init(&argc, &argv);
    int rank, commsize;
    int n = (argc > 1) ? atoi(argv[1]) : 100;
    const char *genname = (argc > 2) ? argv[2] : NULL;
    const char *parname = (argc > 3) ? argv[3] : NULL;
    double t = wtime();
    if (n < 1) {
        fprintf(stderr, "Invalid size of graph");
        exit(EXIT_FAILURE);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    int *g = gen_graph(n);
    if (genname)
        write_matrix(genname, n, g);

    int *gp = shortest_paths_p(n, g);

    printf("n:     %d\n", n);
    printf("Time parallel (sec): %.6f\nS(n): %.6f\n", tp);

    if (sername)
        write_matrix(sername, n, gs);
    if (parname)
        write_matrix(parname, n, gp);

    free(g);
    free(gs);
    free(gp);
    return 0;
}
