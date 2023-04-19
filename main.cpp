#include <iostream>
#include <random>
#include <vector>
#include "TROOT.h"
#include "TGraph.h"
#include "TCanvas.h"

// double const MAX_T = 3.25;
// double const MIN_T = 1.25;
// double const STEP_T = 0.25;
static const double NUM_T = 24;
double TEMPERATURES[] = {1.25, 1.5, 1.65,1.8, 1.9, 1.95, 2.0, 2.05,2.1, 2.15,2.2, 2.25,2.3, 2.35,2.4, 2.45,2.5, 2.55,2.6, 2.7,2.75, 3.0, 3.25, 3.5};
static const double NUM_L = 1;
double L_VALUES[] = {80, 40, 80};
const int NUM_MC = 100000;

void create_graph(std::vector<double> x_values, std::vector<double> y_values, char* title, char* filename) {
    // Plotting average magnetization
    double x[x_values.size()], y[y_values.size()];
    int n = x_values.size();
    for (int k = 0; k < n; k++) {
        x[k] = x_values[k];
        y[k] = y_values[k];
    }

    TCanvas* canvas = new TCanvas(); // Create a new canvas object
    TGraph* g = new TGraph(n, x, y);
    g->SetTitle(title);
    g->SetMarkerStyle(kCircle);
    g->SetMarkerColor(kBlue);
    //g->SetLineColor(kBlue);
    //g->SetLineWidth(2);
    g->Draw("AP");
    canvas->Draw("AP"); // Display the canvas object
    canvas->SaveAs(filename);

    delete canvas;
    delete g;
}

void print_s(int** s, int n) {
    for(int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << s[i][j] << "\t";
        }
        std::cout << "\n";
    }
}

void print_v(double* v, int n) {
    for (int i = 0; i < n; i++) {
        std::cout << v[i] << ", ";
    }
    std::cout << "\n";
}

void print_std_v(std::vector<double> v) {
    for (double i: v) {
        std::cout << i << ", ";
    }
    std::cout << "\n";
}

double average(double* v, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += v[i];
    }
    return sum/n;
}

void unordered_initialize_s(int** s, double L, std::mt19937 GEN, std::uniform_real_distribution<> uniform) {
    for(int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (uniform(GEN) < 0.5) {
                s[i][j] = -1;
            } else {
                s[i][j] = 1;
            }
        }
    }
}

double binder_ratio(double avg_m4, double avg_m2) {
    double br = 1 - avg_m4 / (3 * std::pow(avg_m2,2) );
    return br;
}

double magnetic_susceptibility(double N, double T, double avg_m2, double avg_mgn) {
    double var = avg_m2 - std::pow(avg_mgn, 2);
    double ms = N / T * var; // k considered k = 1
    return ms;
}

double magnetization(int** s, int L, int N) {
    double sum = 0.0;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            sum += s[i][j];
        }
    }
    return sum/N;
}

double* ising_metropolis(int** s, int L, int N,  double T, double* m_values, double* m2_values, double* m4_values,
                         std::mt19937 GEN,
                         std::uniform_real_distribution<> uniform, std::uniform_int_distribution<> select_point_uniform) {
    for (int mc = 0; mc < NUM_MC; mc++) {
        // Montecarlo Step N*N
        for (int r = 0; r < N; r++) {
            int n = select_point_uniform(GEN);
            int m = select_point_uniform(GEN);
            double incrE = 2 * s[n][m] *
                    (s[(n + 1) % L][m] +
                    s[((n - 1) % L + L) % L][m] +
                    s[n][(m + 1) % L] +
                    s[n][((m - 1) % L + L) % L]);
            double p = std::min(1.0, std::exp(-(incrE/T)));

            if (uniform(GEN) < p) {
                s[n][m] = -s[n][m];
            }

        }

        double m =  magnetization(s, L, N);
        m_values[mc] = m;
        m2_values[mc] = std::pow(m, 2);
        m4_values[mc] = std::pow(m, 4);
    }
    return m_values;
}

int main() {

    for (int l = 0; l < NUM_L; l++) {
        // 2D Lattice configuration.
        const int L = L_VALUES[l];
        const int N = L * L;

        // Random variables and generators.
        std::random_device rd;
        std::mt19937 GEN(rd());
        std::uniform_real_distribution<> uniform(0.0, 1.0);
        std::uniform_int_distribution<> select_point_uniform(0, L - 1);

        // Allocation of the 2D Lattice.
        int** s = new int* [L];
        for(int i = 0; i < L; i++) {
            s[i] = new int[L];
        }

        // Declaring variables.
        double* m_values = new double[NUM_MC];
        double* m2_values = new double[NUM_MC];
        double* m4_values = new double[NUM_MC];
        std::vector<double> temperatures;
        std::vector<double> avg_m_values;
        std::vector<double> ms_values;
        std::vector<double> binder_values;

        for (int t = 0; t < NUM_T; t++) {
            const double T = TEMPERATURES[t];
            std::cout << "Computing T = " << T << std::endl;

            unordered_initialize_s(s, L, GEN, uniform);
            ising_metropolis(s, L, N, T, m_values, m2_values, m4_values, GEN, uniform, select_point_uniform);

            // Compute the average magnetization value
            double avg_m = std::abs(average(m_values, NUM_MC));
            // Auxiliary calculations
            double avg_m2 = average(m2_values, NUM_MC);
            double avg_m4 = average(m4_values, NUM_MC);

            // Compute the magnetic susceptibility value
            double ms = magnetic_susceptibility(N, T, avg_m2, avg_m);

            // Compute the binder ratio value
            double br = binder_ratio(avg_m4, avg_m2);
            std::cout << "T = " << T << ", avg_mgn = " << avg_m << ", ms = " << ms << ", br = " << br << std::endl;

            temperatures.push_back(T);
            avg_m_values.push_back(avg_m);
            ms_values.push_back(ms);
            binder_values.push_back(br);
        }

        // Print results.
        std::cout << "L = " << L << ", N = " << N << ", NUM_MC = " << NUM_MC << ", NUM_T = " << NUM_T << std::endl;
        print_std_v(temperatures);
        print_std_v(avg_m_values);
        print_std_v(ms_values);
        print_std_v(binder_values);

        // Create graphs with the results.
        create_graph(temperatures, avg_m_values, "Average Magnetization", "avg_m.png");
        create_graph(temperatures, ms_values, "Magnetic susceptibility", "ms.png");
        create_graph(temperatures, binder_values, "Binder Ratio", "br.png");


        // Deallocation of the 2D Lattice.
        for (int i = 0; i < L; i++) {
            delete[] s[i];
        }
        delete[] s;
        // Deallocation of auxiliary arrays.
        delete[] m4_values;
        delete[] m2_values;
        delete[] m_values;
    }

    return 0;
}

