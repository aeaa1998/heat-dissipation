#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <iostream>
#include <chrono>

// Precisión o diferencia requerida
float ERR = 0.00003f;

// C constant
#define C_CONSTANT 0.5f
// Número de intervalos discretos 10ˆ5
#define DIFUSIVE_VAL 10e-4f
// Length
#define L 20.0f
// Intervals
// int N = 10000;
int N = 1000;

// Delta t, x
float DT,
    DX;
int T;

void resolveInitialValues()
{
    DX = L / N;
    DT = (C_CONSTANT * pow(DX, 2)) / DIFUSIVE_VAL;
    // MIN_EVALUATIONS = N / 2;
    // MIN_EVALUATIONS = 1;
}

int main(int argc, char *argv[])
{

    resolveInitialValues();

    T = DT * N;

    // Temperatures array
    float t0 = 1.0f;
    float tl = 0.0f;
    float tr = 10.0f;
    // x initial temps

    float temperatures[3][N + 1];

    // Set left
    temperatures[2][0] = temperatures[0][0] = tl;
    temperatures[1][0] = temperatures[1][0] = tl;
    temperatures[0][0] = temperatures[2][0] = tl;

    // Set right
    temperatures[2][N] = temperatures[0][N] = tr;
    temperatures[1][N] = temperatures[2][N] = tr;
    temperatures[0][N] = temperatures[1][N] = tr;
    // Set initial temperatures1 excluding the left and right one
    float avg = 0.0f;
    // float evaluation = 0;
    int my_num_threads = 6;
    int temp_pointer = 0;
    temp_pointer++;
    omp_lock_t section_locks_read[3][my_num_threads];
    for (int j = 0; j < 3; j++)
    {
        for (int i = 0; i < 3; i++)
        {
            omp_init_lock(&section_locks_read[j][i]);
            // 0 is always gonna be fine
            if (j != 0 && i != 0)
            {
                omp_set_lock(&section_locks_read[j][i]);
            }
        }
    }

    auto t_start = std::chrono::high_resolution_clock::now();
    int block_size = (N - 2) / my_num_threads;
    float f = -1.0f;
#pragma omp parallel for default(none) shared(temperatures, t0, N)
    for (int i = 1; i < N; i++)
    {
        temperatures[2][i] = t0;
    }
    do
    {
        avg = 0.0f;
        int numberToCheck = 0;
        // i stands for the current x partition in our l bar
#pragma omp parallel for collapse(2) num_threads(my_num_threads) schedule(static, block_size) reduction(+ \
                                                                                                        : avg, numberToCheck) shared(section_locks_read)
        for (int i = 0; i < 3; i++)
        {
            for (int j = 1; j < N; ++j)
            {
                // Only in the first block if it is not last thread
                int thread_id = omp_get_thread_num();
                int l = j - (thread_id * block_size);

                int prev = i == 0 ? 2 : i - 1;
                int next = i == 2 ? 2 : i + 1;

                if (i != 0 && thread_id != 0)
                {
                    // WAIT TO READ AND CLEAN
                    printf("Locked thread %d in pos i %d\n", thread_id, prev);
                    omp_set_lock(&section_locks_read[prev][thread_id]);
                }
                printf("Passed thread %d in pos i %d\n", thread_id, i);
                temperatures[i][j] = temperatures[prev][j] + C_CONSTANT * (temperatures[prev][j - 1] - 2 * temperatures[prev][j] + temperatures[prev][j + 1]);
                // Is the last iteration
                if (i == 2 && temperatures[i][j] != tl && temperatures[i][j] != tr && temperatures[i][j] != t0)
                {
                    avg += temperatures[i][j] - temperatures[prev][j];
                    numberToCheck++;
                }

                if (thread_id != my_num_threads - 1)
                {
                    int thread_to_block = thread_id + 1;
                    // printf("Unlocked thread %d in pos i %d\n", thread_to_block, i);
                    // omp_unset_lock(&section_locks_read[i][thread_to_block]);
                }
            }
        }
        f = temperatures[2][N / 2];
        printf("%f\n", f);
#pragma omp parallel for collapse(2)
        for (int j = 0; j < 3; j++)
        {
            for (int i = 0; i < 3; i++)
            {
                if (j != 0 && i != 0)
                {
                    omp_set_lock(&section_locks_read[j][i]);
                }
            }
        }

    } while (f > (tr / 2) + 0.5f || f < (tr / 2) - 0.5f);

    // for (int i = 0; i < N + 1; i++)
    // {
    //     printf("In position %d value %f\n", i, temperatures[2][i]);
    // }
    printf("In position %d value %f\n", N / 2, temperatures[2][N / 2]);
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    printf("El tiempo que paso fue %f milisegundos", elapsed_time_ms);
    return 0;
}