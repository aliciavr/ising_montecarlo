#include <iostream>
#include <random>
#include <vector>
#include "TROOT.h"

const int N = 100;
double T = 3.0; // T belongs to [0, 5]
double const MAX_T = 4.0;
const int NUM_CHANGES = N * N;
const int NUM_MC = 10000;

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

int** unordered_initialize_s(int** s) {
    for(int i = 0; i < N; i++) {
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

void print_std_v(std::vector<double> v) {
    for (double i: v) {
        std::cout << i << "\t";
    }
    std::cout << "\n";
}

double average(double* v, int N) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += v[i];
    }
    return sum/N;
}

double* ising_metropolis(int** s, double* magnetization_values) {

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
    int** s = new int* [N];
    for(int i = 0; i < N; i++) {
        s[i] = new int[N];
    }
    double* magnetization_values = new double[NUM_MC];
    std::vector<double> avg_mgn_values;
    std::vector<double> temperatures;

    for (double t = 0.5; t < MAX_T; t += 0.25) {
        std::cout << "Computing T = " << t << std::endl;
        unordered_initialize_s(s);
        ising_metropolis(s, magnetization_values);
        //print_s(s, N);
        //print_v(magnetization_values, NUM_MC);
        double avg_magnetization = average(magnetization_values, NUM_MC);
        std::cout << "T = " << t << ", avg_magnetization = " << avg_magnetization << std::endl;
        avg_mgn_values.push_back(avg_magnetization);
        temperatures.push_back(t);

    }

    print_std_v(avg_mgn_values);
    print_std_v(temperatures);


    delete[] magnetization_values;
    for (int i = 0; i < N; i++) {
        delete[] s[i];
    }
    delete[] s;

    return 0;
}
