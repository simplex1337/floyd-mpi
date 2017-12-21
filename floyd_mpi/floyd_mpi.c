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

int min(int a, int b)
{
    return (a < b) ? a : b;
}

void get_chunk(int a, int b, int commsize, int rank, int *lb, int *ub)
{
    int n = b - a + 1;
    int q = n / commsize;
    if (n % commsize)
        q++;
    int r = commsize * q - n;
    /* Compute chunk size for the process */
    int chunk = q;
    if (rank >= commsize - r) chunk = q - 1;
    *lb = a; /* Determine start item for the process */
    if (rank > 0) { /* Count sum of previous chunks */
        if (rank <= commsize - r)
            *lb += q * rank;
        else
            *lb += q * (commsize - r) + (q - 1) * (rank - (commsize - r));
    }
    *ub = *lb + chunk - 1;
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

void copy(int *p_beg, int *p_end, int *rows)
{
    int *p = p_beg;
    for (int i = 0; p < p_end; i++, p++)
        rows[i] = *p;
}

// функция для рассылки строки всем процессам
void raw_distribution(int *proc_raws, int n, int raw_num, int k, int *raw, int commsize, int rank)
{
    int proc_raw_rank = 228; //ранг процесса, которому принадлежит k-я строка
    int proc_raw_num; //номер k-й строки в полосе матрицы

    //нахождение ранга процесса - владельца k-й строки
    /*
    int rest_raws = n;
    int ind = 0;
    int num = n / commsize;
    for (proc_raw_rank = 1; proc_raw_rank < commsize + 1; proc_raw_rank++) {
        if (k < ind + num)
            break;
        rest_raws -= num;
        ind += num;
        num = rest_raws / (commsize - proc_raw_rank);
    }
    */
    int lb, ub;
    for (int i = 0; i < commsize; i++) {
        get_chunk(0, n - 1, commsize, i, &lb, &ub);
        if (k >= lb && k <= ub) {
            proc_raw_rank = i;
            break;
        }
    }

    //копирование строки  в массив raw
    if (proc_raw_rank == rank) {
        int need_row = k - lb;
        copy(&proc_raws[need_row * n], &proc_raws[(need_row + 1) * n], raw);
    }
    //распределение k-й строки между процессами
    MPI_Bcast(raw, n, MPI_INT, proc_raw_rank, MPI_COMM_WORLD);
}

void floyd_mpi(int *proc_raws, int n, int raw_num, int commsize, int rank)
{
    int *raw = (int *) malloc(sizeof(*raw) * n);
    int t1, t2;

    for (int k = 0; k < n; k++) {
        //распределение k-й строки среди процессов
        raw_distribution(proc_raws, n, raw_num, k, raw, commsize, rank);
        // обновление  элементов матрицы смежности
        for (int i = 0; i < raw_num; i++) {
            for (int j = 0; j < n; j++) {
                if( (proc_raws[i * n + k] != -1) && (raw[j] != -1)) {
                    t1 = proc_raws[i * n + j];
                    t2 = proc_raws[i * n + k] + raw[j];
                    proc_raws[i * n + j] = min(t1, t2);
                }
            }
        }
    }
    free(raw);
}

//распределение данных между процессами
void data_distribution(int *g, int *proc_raws, int n, int raw_num)
{

}

void result_collection(int *g, int *proc_raws, int raw_num)
{

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
    int *g = gen_graph(n); //инициализация
    int rank, commsize;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    //if (genname && !rank)
        //write_matrix(genname, n, g);

    MPI_Bcast(g, n * n, MPI_INT, 0, MPI_COMM_WORLD);

    /*if (rank == 0) {
        for (int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++)
                printf("%d ", g[i * n + j]);
            printf("\n");
        }
    }*/

    int lb, ub;
    get_chunk(0, n - 1, commsize, rank, &lb, &ub);
    int raw_num = ub - lb + 1; //чисто строк для текущего процесса
    //строки матрицы смежности текущего процесса
    int *proc_raws = malloc(sizeof(*proc_raws) * raw_num * n);
    for (int i = 0; i < raw_num; i++)
        for (int j = 0; j < n; j++)
            proc_raws[i * n + j] = (g[(i + lb) * n + j] == 0) ? 10000 : g[(i + lb) * n + j];

    // int *gp = (int *) calloc(n * n, sizeof(int));
    // memcpy(gp, g, n * n * sizeof(int)); //матрица хранится в каждом процессе
    // infinitize(n, proc_raws);


    floyd_mpi(proc_raws, n, raw_num, commsize, rank);

    // deinfinitize(n, proc_raws);

    if (!rank) {
        for (int i = 0; i < n * n; i++)
            g[i] = 0;
        int *recvcnt = malloc(sizeof(*recvcnt) * commsize);
        int *displs = malloc(sizeof(*displs) * commsize);
        for (int i = 0; i < commsize; i++) {
            int lb, ub;
            get_chunk(0, n - 1, commsize, i, &lb, &ub);
            recvcnt[i] = (ub - lb + 1) * n;
            displs[i] = (i > 0) ? displs[i - 1] + recvcnt[i - 1] : 0;
            //printf("cnt = %d, displc = %d\n", recvcnt[i], displs[i]);
        }
        MPI_Gatherv(proc_raws, raw_num * n, MPI_INT, g, recvcnt, displs, MPI_INT, 0, MPI_COMM_WORLD);

    } else
        MPI_Gatherv(proc_raws, raw_num * n, MPI_INT, NULL, NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < n * n; i += n + 1)
       g[i] = 0;
    t = wtime() - t;
    double tmax;
    MPI_Reduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!rank) {
        printf("n:     %d\n", n);
        printf("Time for %d proc: %.6f\n", commsize, tmax);
    }
    /*if (rank == 0) {
        printf("----\n");
        for (int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++)
                printf("%d ", g[i * n + j]);
            printf("\n");
        }
    }
    //if (parname)
      //  write_matrix(parname, n, gp);
    */
    free(g);
    free(proc_raws);
    MPI_Finalize();
    return 0;
}
