#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>

// Precisión o diferencia requerida
double ERR = 0.00001;
// C constant
#define C_CONSTANT 0.5
// Número de intervalos discretos 10ˆ5
#define DIFUSIVE_VAL 10e-4f
// Length
#define L 10.0f
// Intervals
int N = 2000;

// Delta t, x
double DT, DX;
int T;

double maxFromArray(double *data)
{
    double max = 0;
    int i;
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
    // MIN_EVALUATIONS = N / 2;
    // MIN_EVALUATIONS = 1;
}

int main(int argc, char *argv[])
{

    resolveInitialValues();

    T = DT * N;

    // Temperatures array
    double t0 = 60.0f;
    double tl = 100.0f;
    double tr = 40.0f;
    // x initial temps

    double temperatures1[N + 1];
    double temperatures2[N + 1];
    double differences[N + 1];
    // Set left
    temperatures1[0] = temperatures2[0] = tl;
    // Set right
    temperatures1[N] = temperatures2[N] = tr;
    // Set initial temperatures1 excluding the left and right one
    // double evaluation = 0;

    for (int i = 1; i < N; i++)
    {
        temperatures1[i] = t0;
    }
    double start_time = omp_get_wtime();
    double max = 1;
    while (max > ERR)
    {
        // i stands for the current x partition in our l bar
        for (int i = 1; i < N; i++)
        {
            temperatures2[i] = temperatures1[i] + C_CONSTANT * (temperatures1[i - 1] - 2 * temperatures1[i] + temperatures1[i + 1]);
            differences[i] = fabs(temperatures2[i] - temperatures1[i]);
        }

        memcpy(temperatures1, temperatures2, sizeof(temperatures1));

        double new_max = maxFromArray(differences);

        if (new_max < ERR)
        {
            max = new_max;
        }
    }

    double elapsed_time_ms = omp_get_wtime() - start_time;
    printf("El tiempo que paso fue en el cálculo de temperaturas %f segundos \n", elapsed_time_ms);

    printf("In position %d value %f\n", 0, temperatures1[0]);
    printf("In position %d value %f\n", N / 5, temperatures1[N / 5]);
    printf("In position %d value %f\n", (N / 5) * 2, temperatures1[(N / 5) * 2]);
    printf("In position %d value %f\n", (N / 5) * 3, temperatures1[(N / 5) * 3]);
    printf("In position %d value %f\n", (N / 5) * 4, temperatures1[(N / 5) * 4]);
    printf("In position %d value %f\n", (N / 5) * 5, temperatures1[(N / 5) * 5]);
    printf("In position %d value %f\n", N / 2, temperatures1[N / 2]);

    return 0;
}