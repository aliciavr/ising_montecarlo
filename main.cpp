#include <iostream>
#include <random>
#include <vector>
#include "TROOT.h"
#include "TGraph.h"
#include "TCanvas.h"

const int N = 20;
double const MAX_T = 3.25;
double const MIN_T = 1.25;
double const STEP_T = 0.15;
const int NUM_CHANGES = N * N;
const int NUM_MC = 100000;

// Random variables
std::random_device rd;
std::mt19937 GEN(rd());
std::uniform_real_distribution<> uniform(0.0, 1.0);
std::uniform_int_distribution<> select_point_uniform(0, N - 1);

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

void unordered_initialize_s(int** s) {
    for(int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (uniform(GEN) < 0.5) {
                s[i][j] = -1;
            } else {
                s[i][j] = 1;
            }
        }
    }
}

double magnetic_susceptibility(double T, double* magnetization2_values, double avg_mgn) {
    double var = average(magnetization2_values, NUM_MC) - std::pow(avg_mgn, 2);
    double ms = N / T * var; // k considered k = 1
    return ms;
}

double magnetization(int** s) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            sum += s[i][j];
        }
    }
    return sum/(N*N);
}

double* ising_metropolis(int** s, double T, double* magnetization_values, double* magnetization2_values) {

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

        double m =  magnetization(s);
        magnetization_values[mc] = m;
        magnetization2_values[mc] = m * m;
    }
    return magnetization_values;
}

int main() {
    int** s = new int* [N];
    for(int i = 0; i < N; i++) {
        s[i] = new int[N];
    }
    double* magnetization_values = new double[NUM_MC];
    double* magnetization2_values = new double[NUM_MC];
    std::vector<double> temperatures;
    std::vector<double> avg_mgn_values;
    std::vector<double> ms_values;
    for (double T = MIN_T; T < MAX_T; T += STEP_T) {
        std::cout << "Computing T = " << T << std::endl;
        unordered_initialize_s(s);
        ising_metropolis(s, T, magnetization_values, magnetization2_values);
        //print_s(s, N);
        //print_v(magnetization_values, NUM_MC);
        // Store the temperature value
        temperatures.push_back(T);
        // Store the average magnetization value
        double avg_mgn = std::abs(average(magnetization_values, NUM_MC));
        avg_mgn_values.push_back(avg_mgn);
        // Store the magnetic susceptibility value
        double ms = magnetic_susceptibility(T, magnetization2_values, avg_mgn);
        ms_values.push_back(ms);
        std::cout << "T = " << T << ", avg_mgn = " << avg_mgn << ", ms = " << ms << std::endl;

    }

    print_std_v(temperatures);
    print_std_v(avg_mgn_values);
    print_std_v(ms_values);

    // Plotting average magnetization
    double x[temperatures.size()], y[avg_mgn_values.size()];
    int n = temperatures.size();
    for (int k = 0; k < n; k++) {
        x[k] = temperatures[k];
        y[k] = avg_mgn_values[k];
    }

    TCanvas* canvas = new TCanvas(); // Create a new canvas object
    TGraph* g = new TGraph(n, x, y);
    g->SetTitle("Magnetizaci\'on promedio en funci\'on de la temperatura; Temperaturas; Magnetizaci\'on promedio");
    g->SetMarkerStyle(kCircle);
    g->SetMarkerColor(kBlue);
    //g->SetLineColor(kBlue);
    //g->SetLineWidth(2);
    g->Draw("AP");
    canvas->Draw("AP"); // Display the canvas object
    canvas->SaveAs("graph.png");

    // Plotting magnetic susceptibility


    delete[] magnetization2_values;
    delete[] magnetization_values;
    for (int i = 0; i < N; i++) {
        delete[] s[i];
    }
    delete[] s;

    return 0;
}
