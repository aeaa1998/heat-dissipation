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
int N = 10000;

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
    float t0 = 25.f;
    float tl = 0.0f;
    float tr = 100.0f;
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
    auto t_start = std::chrono::high_resolution_clock::now();
    int my_num_threads = 6;
    int temp_pointer = 0;
    temp_pointer++;
    omp_lock_t section_locks_read[3][my_num_threads];
    omp_lock_t section_locks_write[3][my_num_threads];
    int block_size = (N - 2) / my_num_threads;

#pragma omp parallel for default(none) shared(temperatures, t0, N)
    for (int i = 1; i < N; i++)
    {
        temperatures[2][i] = t0;
    }
    do
    {
        avg = 0.0f;
        int numberToCheck = 0;
        int counter = 0;
        // i stands for the current x partition in our l bar
#pragma omp parallel num_threads(my_num_threads) reduction(+ \
                                                           : avg, numberToCheck, counter)
        for (int i = 0; i < 3; i++)
        {

            for (int d = 0; d < block_size; ++d)
            {
                // Only in the first block if it is not last thread
                int thread_id = omp_get_thread_num();
                int j = ((block_size + 1) * thread_id) + d;
                if (j == 1 && thread_id != my_num_threads - 1)
                {
                    int thread_to_block = thread_id + 1;
                    omp_set_lock(&section_locks_read[i][thread_to_block]);
                }

                if (i != 0 && thread_id != 0)
                {
                    // WAIT TO READ AND CLEAN
                    omp_set_lock(&section_locks_read[i][thread_id]);
                    omp_unset_lock(&section_locks_read[i][thread_id]);
                }

                int prev = i == 0 ? 2 : i - 1;

                temperatures[i][j] = temperatures[prev][j] + C_CONSTANT * (temperatures[prev][j - 1] - 2 * temperatures[prev][j] + temperatures[prev][j + 1]);
                // Is the last iteration
                if (j == N - 1)
                {
                    printf("%d %d\n", thread_id, j);
                }

                if (i == 2)
                {
                    counter++;
                    avg += temperatures[i][j] - temperatures[prev][j];
                    numberToCheck++;
                }

                if (j == 1 && thread_id != my_num_threads - 1)
                {
                    int thread_to_block = thread_id + 1;
                    omp_unset_lock(&section_locks_read[i][thread_to_block]);
                }
            }
        }

        if (numberToCheck > 1)
        {
            avg = avg / numberToCheck;
        }
        else
        {
            avg = 1.0f;
        }
        // printf("Avg %d\n", avg);
    } while (avg > ERR);

    for (int i = 0; i < N + 1; i++)
    {
        printf("In position %d value %f\n", i, temperatures[2][i]);
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    printf("B size %d\n", block_size);
    printf("El tiempo que paso fue %f milisegundos", elapsed_time_ms);
    return 0;
}