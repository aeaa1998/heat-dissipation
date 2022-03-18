#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// Precisión o diferencia requerida
float ERR = 0.000005f;
// C constant
#define C_CONSTANT 0.5f
// Número de intervalos discretos 10ˆ5
#define DIFUSIVE_VAL 10e-4f
// Length
#define L 10.0f
// Intervals
int N = 1000;

// Delta t, x
float DT, DX;
int T;

float maxFromArray(float *data)
{
    float max = 0;
    int i;
#pragma omp for private(i)
    for (i = 1; i < N; i++)
    {
        if (max < data[i])
        {
            max = data[i];
        }
    }
    return max;
}

void resolveInitialValues()
{
    DX = L / N;
    DT = (C_CONSTANT * pow(DX, 2)) / DIFUSIVE_VAL;
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
    float *temperatures1;
    float *temperatures2;
    float *differences;
    temperatures1 = new float[N + 1];
    temperatures2 = new float[N + 1];
    differences = new float[N + 1];
    // Set left
    temperatures1[0] = temperatures2[0] = tl;
    // Set right
    temperatures1[N] = temperatures2[N] = tr;
    // Set initial temperatures1 excluding the left and right one
    // float evaluation = 0;
    for (int i = 1; i < N; i++)
    {
        temperatures1[i] = t0;
    }
    auto t_start = std::chrono::high_resolution_clock::now();
    float max = 1, new_max;
    int j, k, m;

    // We create the parallel threads in the while block so they live all the time the while block lives
    // temperatures are shared as N, tr and res
#pragma omp parallel shared(temperatures1, temperatures2, N, tr, differences, ERR) private(j, k, m, new_max, max)
    while (max > ERR)
    {

#pragma omp parallel for
        for (j = 1; j < N; j++)
        {
            temperatures2[j] = temperatures1[j] + C_CONSTANT * (temperatures1[j - 1] - 2 * temperatures1[j] + temperatures1[j + 1]);
            differences[j] = fabs(temperatures2[j] - temperatures1[j]);
        }

#pragma omp parallel for
        for (k = 1; k < N; k++)
        {
            temperatures1[k] = temperatures2[k];
        }
        new_max = 0;
#pragma omp parallel for
        for (m = 1; m < N; m++)
        {
            if (new_max < differences[m])
            {
                new_max = differences[m];
            }
        }
        if (new_max < ERR)
        {
            // max = new_max;
            // #pragma omp flush(max)
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    printf("El tiempo que paso fue en el cálculo de temperaturas %f milisegundos \n", elapsed_time_ms);
    printf("In position %d value %f\n", 0, temperatures1[0]);
    printf("In position %d value %f\n", N / 5, temperatures1[N / 5]);
    printf("In position %d value %f\n", (N / 5) * 2, temperatures1[(N / 5) * 2]);
    printf("In position %d value %f\n", (N / 5) * 3, temperatures1[(N / 5) * 3]);
    printf("In position %d value %f\n", (N / 5) * 4, temperatures1[(N / 5) * 4]);
    printf("In position %d value %f\n", (N / 5) * 5, temperatures1[(N / 5) * 5]);
    printf("In position %d value %f\n", N / 2, temperatures1[N / 2]);
    return 0;
}