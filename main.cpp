#include <iostream>
#include <random>
#include "TROOT.h"
int** unordered_initialize_s(const int N, std::mt19937 gen, std::uniform_real_distribution<> uniform) {
    int** s = new int* [N];
    for(int i = 0; i < N; i++) {
        s[i] = new int[N];
        for (int j = 0; j < N; j++) {
            if (uniform(gen) < 0.5) {
                s[i][j] = -1;
            } else {
                s[i][j] = 1;
            }
        }
    }
    return s;
}

void print_s(int** s, const int N) {
    for(int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << s[i][j] << "\t";
        }
        std::cout << "\n";
    }
}

void ising_metropolis(int** s, const int N, const double T, std::mt19937 gen, std::uniform_real_distribution<> uniform) {
    //int** s = unordered_initialize_s(N, gen, uniform);

}
int main() {
    // Parameters
    const int N = 10;
    double T = 2.0; // T belongs to [0, 5]

    // Random variables
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    std::uniform_int_distribution<> select_point_uniform(0, N - 1);

    // Unordered initialization of the Ising Matrix.
    int changed[N][N] = {};
    int** s = unordered_initialize_s(N, gen, uniform);
    print_s(s, N);


    /*
    int n = select_point_uniform(gen);
    int m = select_point_uniform(gen);

    std::cout << "n = " << n << " m = " << m << std::endl;

    double incrE = 2 * s[n][m] * (s[(n + 1) % N][m] + s[(n - 1) % N][m] + s[n][(m + 1) % N] + s[n][(m - 1) % N]);
    double p = std::min(1.0, std::exp(-(incrE/T)));

    std::cout << "p = " << p << std::endl;

    if (uniform(gen) < p) {
        s[n][m] = -s[n][m];
        changed[n][m] = 1;
    }
    */
    /*
    for (int r = 0; r < N * N; r++) {

    }
     */
    for (int i = 0; i < N; i++) {
        delete[] s[i];
    }
    delete[] s;

    return 0;
}
