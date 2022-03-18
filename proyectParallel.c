#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <chrono>
#include <omp.h>
#include <string.h>

// Precisión o diferencia requerida
#define ERR 0.00001

// C constant
#define C_CONSTANT 0.5
// Número de intervalos discretos 10ˆ5
#define DIFUSIVE_VAL 10e-4
// Length
#define L 10.0
// Intervals
#define N 5000
// Delta t, x
double DT, DX;
int T;

double maxFromArray(double data[], int start, int end)
{
    double inner_max = 0;
    int i;
    // #pragma omp parallel for shared(data, end, start) private(i) reduction(max \
//                                                                        : inner_max)
    for (i = start; i < end; i++)
    {
        inner_max = inner_max > data[i] ? inner_max : data[i];
    }
    return inner_max;
}

void resolveInitialValues()
{
    DX = L / N;
    DT = (C_CONSTANT * pow(DX, 2)) / DIFUSIVE_VAL;
}

void sendTask(int start, int end, int num_threads, double temperatures2[], double temperatures1[], double diff[])
{
    double max = 0;

    // #pragma omp task
    for (int j = start; j < end; j++)
    {
        // if (j == end - 1)
        // {
        //     printf("here %d \n", j);
        // }
        temperatures2[j] = temperatures1[j] + C_CONSTANT * (temperatures1[j - 1] - 2 * temperatures1[j] + temperatures1[j + 1]);
        diff[j] = fabs(temperatures2[j] - temperatures1[j]);
    }
}

int main(int argc, char *argv[])
{

    resolveInitialValues();

    T = DT * N;
    double t0 = 60.0;
    double tl = 100.0;
    double tr = 40.0;
    int num_threads = 4;
    if (argc > 4)
    {
        num_threads = strtol(argv[1], NULL, 10);
        t0 = strtod(argv[2], NULL);
        tl = strtod(argv[3], NULL);
        tr = strtod(argv[4], NULL);
    }

    // x initial temps
    double temperatures1[N];
    double temperatures2[N];

    // Set left
    temperatures1[0] = tl;
    temperatures2[0] = tl;
    // Set right
    temperatures1[N - 1] = tr;
    temperatures2[N - 1] = tr;

    for (int i = 1; i < N - 1; i++)
    {
        temperatures1[i] = t0;
    }

    double max = 1;
    double new_max;
    int m;
    int iterations = (N - 2);
    int block = ceil((float)(iterations) / (float)(num_threads));

    // if ((block * num_threads) < iterations)
    // {
    // block -= 1;
    // }
    double differences[N];
    double differences_pt[num_threads];
    // printf("B size: %d", block);
    double start_time = omp_get_wtime();
    omp_set_num_threads(num_threads);
    omp_set_dynamic(0);
    int end_l = (iterations - ((num_threads - 1) * block)) % iterations;
#pragma omp parallel shared(temperatures2, temperatures1, differences, max, new_max, block, num_threads, m, differences_pt, iterations, end_l)
#pragma omp single
    {
        while (max > ERR)
        {
            // Al threads create just one part of the block
            for (m = 0; m < num_threads - 1; m++)
            {
#pragma omp task shared(temperatures2, temperatures1, differences, block, num_threads, m, differences_pt, iterations, end_l)
                {
                    int start = m * block + 1;
                    int end = start + block;
                    // printf("Task %d started with %d ended with %d with b %d\n", m, start, end, block);
                    sendTask(start, end, num_threads, temperatures2, temperatures1, differences);
                    differences_pt[m] = maxFromArray(differences, start, end);
                }
                // printf("T%d done\n", tid);
            }
#pragma omp task shared(temperatures2, temperatures1, differences, block, num_threads, differences_pt, iterations, end_l)
            {

                int start = (num_threads - 1) * block + 1;
                int end = start + end_l;
                // printf("Task %d started with %d ended with %d with b %d\n", m, start, end, block);
                sendTask(start, end, num_threads, temperatures2, temperatures1, differences);
                differences_pt[m] = maxFromArray(differences, start, end);
            }
            // printf("T%d done\n", tid);

#pragma omp taskwait

            // printf("T%d started\n", 1);
            // TODO check how to make this parallell
            double new_max = maxFromArray(differences_pt, 0, num_threads);
            // double new_max = maxFromArray(differences, 1, iterations);
            // printf("all done\n");
            memcpy(temperatures1, temperatures2, sizeof(temperatures1));
            // printf("Nuevo mayor: %lf\n", new_max);
            if (new_max < ERR)
            {
                max = new_max;
            }
        }
        // printf("T%d unlocking\n", 1);
    }

    double elapsed_time_ms = omp_get_wtime() - start_time;
    printf("El tiempo que paso fue en el cálculo de temperaturas %f segundos \n", elapsed_time_ms);

    printf("In position %d value %f\n", 0, temperatures1[0]);
    printf("In position %d value %f\n", N / 5, temperatures1[N / 5]);
    printf("In position %d value %f\n", (N / 5) * 2, temperatures1[(N / 5) * 2]);
    printf("In position %d value %f\n", (N / 5) * 3, temperatures1[(N / 5) * 3]);
    printf("In position %d value %f\n", (N / 5) * 4, temperatures1[(N / 5) * 4]);
    printf("In position %d value %f\n", ((N / 5) * 5) - 1, temperatures1[((N / 5) * 5) - 1]);
    printf("In position %d value %f\n", N / 2, temperatures1[N / 2]);
    return 0;
}