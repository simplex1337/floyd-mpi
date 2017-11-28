#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <sys/time.h>

double tp;

double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double) t.tv_sec + (double) t.tv_usec * 1E-6;
}

void raw_distribution(*proc_raws, int size, int raw_num, int k, int *p_raw)
{

}

void floyd_mpi(int *proc_raws, int size, int raw_num)
{

}

void data_distribution(int *g, *proc_raws, int size, int raw_num)
{

}

void result_collection(int *g, int proc_raws, int raw_num)
{

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
    const char *parname = (argc > 3) ? argv[3] : NULL;
    double t = wtime();
    if (n < 1) {
        fprintf(stderr, "Invalid size of graph");
        exit(EXIT_FAILURE);
    }
    int *g = gen_graph(n);
    if (genname)
        write_matrix(genname, n, g);
    int rank, commsize;
    int *proc_raws; //строки матрицы смежности текущего процесса
    int raw_num; //чисто строк для текущего процесса
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    int *gp = (int *) calloc(n * n, sizeof(int));
    memcpy(gp, g, n * n * sizeof(int));
    infinitize(n, gp);

    //my code here

    for (int i = 0; i < n * n; i += n + 1)
        gp[i] = 0;
    deinfinitize(n, gp);

    // printf("n:     %d\n", n);
    // printf("Time parallel (sec): %.6f\nS(n): %.6f\n", tp);

    if (parname)
        write_matrix(parname, n, gp);

    free(g);
    free(gp);
    MPI_Finalize();
    return 0;
}
