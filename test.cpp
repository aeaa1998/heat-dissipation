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

#pragma omp parallel for collapse(2) num_threads(6) schedule(static, 4)
    for (int i = 0; i < 3; i++)
    {
        for (int d = 0; d < 4 * 6; ++d)
        {
            printf("El thread %d esta corriendo el valor i %d y d %d \n", omp_get_thread_num(), i, d);
        }
    }

    return 0;
}