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
    float temperatures1[N + 1];
    float temperatures2[N + 1];

    // Set left
    temperatures1[0] = temperatures2[0] = tl;
    // Set right
    temperatures1[N] = temperatures2[N] = tr;
    // Set initial temperatures1 excluding the left and right one
    float avg = 0.0f;
    // float evaluation = 0;
    auto t_start = std::chrono::high_resolution_clock::now();
    int my_num_threads = 12;

    omp_lock_t section_locks_read[my_num_threads];
    omp_lock_t section_locks_write[my_num_threads];

#pragma omp parallel for default(none) shared(temperatures1, t0, N)
    for (int i = 1; i < N; i++)
    {
        temperatures1[i] = t0;
    }

    do
    {
        avg = 0.0f;
        int numberToCheck = 0;
        // i stands for the current x partition in our l bar

#pragma omp parallel for reduction(+ \
                                   : avg, numberToCheck)
        for (int i = 1; i < N; i++)
        {
            temperatures2[i] = temperatures1[i] + C_CONSTANT * (temperatures1[i - 1] - 2 * temperatures1[i] + temperatures1[i + 1]);

            if (temperatures2[i] != tl && temperatures2[i] != tr && temperatures2[i] != t0)
            {

                avg += temperatures2[i] - temperatures1[i];
                numberToCheck++;
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

#pragma omp parallel for default(none) shared(temperatures1, temperatures2, N)
        for (int i = 1; i < N; i++)
        {
            temperatures1[i] = temperatures2[i];
        }
    } while (avg > ERR);

    for (int i = 0; i < N + 1; i++)
    {
        printf("In position %d value %f\n", i, temperatures1[i]);
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    printf("El tiempo que paso fue %f milisegundos", elapsed_time_ms);
    return 0;
}