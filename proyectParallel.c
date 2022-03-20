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
    // Iterate the array to get current max
    for (i = start; i < end; i++)
    {
        inner_max = inner_max > data[i] ? inner_max : data[i];
    }
    return inner_max;
}

// Resolve current initial values
void resolveInitialValues()
{
    DX = L / N;
    DT = (C_CONSTANT * pow(DX, 2)) / DIFUSIVE_VAL;
}

void sendTask(int start, int end, int num_threads, double temperatures2[], double temperatures1[], double diff[])
{
    double max = 0;

    for (int j = start; j < end; j++)
    {
        // We get the new temperature
        temperatures2[j] = temperatures1[j] + C_CONSTANT * (temperatures1[j - 1] - 2 * temperatures1[j] + temperatures1[j + 1]);
        // And saave the difference
        diff[j] = fabs(temperatures2[j] - temperatures1[j]);
    }
}

int main(int argc, char *argv[])
{

    resolveInitialValues();
    // Resolve the time
    T = DT * N;
    // The initial temperature with a default value
    double t0 = 60.0;
    // The temperature at the left
    double tl = 100.0;
    // Temperatura de la barra a la derecha
    double tr = 40.0;
    // Núm
    int num_threads = 4;
    // if we have the correct number of arguments we set the values
    if (argc > 4)
    {
        num_threads = strtol(argv[1], NULL, 10);
        t0 = strtod(argv[2], NULL);
        tl = strtod(argv[3], NULL);
        tr = strtod(argv[4], NULL);
    }

    // x initial temps
    // Array to store the first
    double temperatures1[N];
    // Array to store the calculated temperatures
    double temperatures2[N];

    // Set left temp on the bar on both arrays because we wont touch it
    temperatures1[0] = tl;
    temperatures2[0] = tl;
    // Set right temp on the bar on both arrays because we wont touch it later
    temperatures1[N - 1] = tr;
    temperatures2[N - 1] = tr;

    // Set the initial temperatures on the first array
    for (int i = 1; i < N - 1; i++)
    {
        temperatures1[i] = t0;
    }

    // Current max value
    double max = 1;
    // Set the new max as a double without a value
    double new_max;
    // Variable m for the iteration of threads
    int m;
    // Number of iterations is N - 2 because we start at 1 and finish at N-2
    int iterations = (N - 2);
    // The block size
    int block = ceil((float)(iterations) / (float)(num_threads));

    // The differences equivalent to the array of temperatures
    double differences[N];
    // The diferences each task will give us
    double differences_pt[num_threads];
    // Get start time
    double start_time = omp_get_wtime();
    // Set the number of threads
    omp_set_num_threads(num_threads);
    // Force it to dont be dynamic
    omp_set_dynamic(0);
    // Get the last iteration for the last thread so it does not have to calculate each time
    int end_l = (iterations - ((num_threads - 1) * block)) % iterations;
    // Create the threads ouside the while
#pragma omp parallel shared(temperatures2, temperatures1, differences, max, new_max, block, num_threads, m, differences_pt, iterations, end_l)
// Run the while in a single thread
#pragma omp single
    {
        while (max > ERR)
        {
            // Al threads create just one part of the block
            // We run all threads excep the last one here
            for (m = 0; m < num_threads - 1; m++)
            {
#pragma omp task shared(temperatures2, temperatures1, differences, block, num_threads, m, differences_pt, iterations, end_l)
                {
                    // Get the start of the thread
                    int start = m * block + 1;
                    // Get the end of the thread
                    int end = start + block;
                    // Calculate the new temps and diffs
                    sendTask(start, end, num_threads, temperatures2, temperatures1, differences);
                    // Get the max diff in this thread
                    differences_pt[m] = maxFromArray(differences, start, end);
                }
            }
            // Execute one last thread but the computation on the first and last values is much cleaner and faster
#pragma omp task shared(temperatures2, temperatures1, differences, block, num_threads, differences_pt, iterations, end_l)
            {
                // Get the start of the thread
                int start = (num_threads - 1) * block + 1;
                // Get the end of the thread
                int end = start + end_l;
                // Calculate the new temps and diffs
                sendTask(start, end, num_threads, temperatures2, temperatures1, differences);
                // Get the max diff in this last thread
                differences_pt[m] = maxFromArray(differences, start, end);
            }

#pragma omp taskwait
            // Explicit barrier telling we need to wait all the tasks

            // Get the real new max
            double new_max = maxFromArray(differences_pt, 0, num_threads);
            // Copy the new temps
            memcpy(temperatures1, temperatures2, sizeof(temperatures1));
            // If the max is bigger than the ERR we set it
            if (new_max < ERR)
            {
                max = new_max;
            }
        }
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