#include <iostream>
#include <random>
#include "TROOT.h"

const int N = 10;
double T = 3.0; // T belongs to [0, 5]
const int NUM_CHANGES = N * N;
const int NUM_MC = 1000;

// Random variables
std::random_device rd;
std::mt19937 GEN(rd());
std::uniform_real_distribution<> uniform(0.0, 1.0);
std::uniform_int_distribution<> select_point_uniform(0, N - 1);

double magnetization(int** s) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            sum += s[i][j];
        }
    }
    return sum/(N*N);
}

int** unordered_initialize_s() {
    int** s = new int* [N];
    for(int i = 0; i < N; i++) {
        s[i] = new int[N];
        for (int j = 0; j < N; j++) {
            if (uniform(GEN) < 0.5) {
                s[i][j] = -1;
            } else {
                s[i][j] = 1;
            }
        }
    }
    return s;
}

void print_s(int** s, int N) {
    for(int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << s[i][j] << "\t";
        }
        std::cout << "\n";
    }
}

void print_v(double* v, int N) {
    for (int i = 0; i < N; i++) {
        std::cout << v[i] << "\t";
    }
    std::cout << "\n";
}

double* ising_metropolis(int** s) {
    double* magnetization_values = new double[NUM_MC];

    for (int mc = 0; mc < NUM_MC; mc++) {
        // Montecarlo Step N*N
        for (int r = 0; r < NUM_CHANGES; r++) {
            int n = select_point_uniform(GEN);
            int m = select_point_uniform(GEN);
            double incrE = 2 * s[n][m] *
                    (s[(n + 1) % N][m] +
                    s[((n - 1) % N + N) % N][m] +
                    s[n][(m + 1) % N] +
                    s[n][((m - 1) % N + N) % N]);
            double p = std::min(1.0, std::exp(-(incrE/T)));

            if (uniform(GEN) < p) {
                s[n][m] = -s[n][m];
            }

        }

        magnetization_values[mc] = magnetization(s);
    }
    return magnetization_values;
}
int main() {

    int** s = unordered_initialize_s();
    double* magnetization_values = ising_metropolis(s);
    print_s(s, N);
    print_v(magnetization_values, NUM_MC);

    delete[] magnetization_values;
    for (int i = 0; i < N; i++) {
        delete[] s[i];
    }
    delete[] s;


    /*
    double* test = new double[4];
    for (int i = 0; i < 4; i++) {
        test[i] = (i + 1) * 10;
    }

    int A = 4;
    print_v(test, A);
    int n = (4 % A + A) % A;
    std::cout << "n = " << n << " test value = " << test[n];

    delete[] test;
     */
    return 0;
}
