#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <chrono>

// Precisión o diferencia requerida
float ERR = 0.0000005f;
// float ERR = 0.05f;

// C constant
#define C_CONSTANT 0.5f
// Número de intervalos discretos 10ˆ5
#define DIFUSIVE_VAL 10e-4f
// Length
#define L 10.0f
// Intervals
int N = 2000;

// Delta t, x
float DT, DX;
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
    float t0 = 15.0f;
    float tl = 20.0f;
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
    float res = -1.0f;
    for (int i = 1; i < N; i++)
    {
        temperatures1[i] = t0;
    }

    do
    {
        avg = 0.0f;
        int numberToCheck = 0;
        // i stands for the current x partition in our l bar
        for (int i = 1; i < N; i++)
        {
            float res = temperatures1[i - 1] - 2 * temperatures1[i] + temperatures1[i + 1];
            temperatures2[i] = temperatures1[i] + C_CONSTANT * res;

            if (temperatures2[i] != tl && temperatures2[i] != tr && temperatures2[i] != t0 && temperatures2[i] - temperatures1[i] != 0)
            {

                avg += temperatures2[i] - temperatures1[i];
                numberToCheck++;
            }
        }

        for (int i = 1; i < N; i++)
        {
            temperatures1[i] = temperatures2[i];
        }
        res = temperatures1[N / 2];
    } while (res > (tr / 2) + 0.1f || res < (tr / 2) - 0.1f);

    for (int i = 0; i < N + 1; i++)
    {
        printf("In position %d value %f\n", i, temperatures1[i]);
    }
    printf("In position %d value %f\n", N / 2, temperatures1[N / 2]);
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    printf("El tiempo que paso fue %f milisegundos", elapsed_time_ms);
    return 0;
}