#include <iostream>
#include <random>
#include <vector>
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"

static const double Tc = 2 / std::log(1 + std::sqrt(2));
static const double NUM_T = 27;
double TEMPERATURES[] = {1.0, 1.1,1.2, 1.3,1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.15,2.2, 2.25,2.3, 2.35, 2.4, 2.45,2.5,2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4};
static const double NUM_L = 4;
double L_VALUES[] = { 10, 20, 40,  80, 100};
const int NUM_MC = 250000;

void create_graph_err(std::vector<double> x_values, std::vector<double> y_values, std::vector<double> err_values,
                      char* title, char* x_label, char* y_label, char* filename) {
    // Plotting average magnetization
    double x[x_values.size()], y[y_values.size()], ey[err_values.size()];
    int n = x_values.size();
    for (int k = 0; k < n; k++) {
        x[k] = x_values[k];
        y[k] = y_values[k];
        ey[k] = err_values[k];
    }


    TCanvas* canvas = new TCanvas(); // Create a new canvas object
    TGraphErrors* g = new TGraphErrors(n, x, y, nullptr, ey);
    g->SetTitle(title);
    g->GetXaxis()->SetTitle(x_label);
    g->GetYaxis()->SetTitle(y_label);
    g->SetMarkerStyle(kCircle);
    g->SetMarkerSize(0.75);
    g->SetMarkerColor(kBlue);
    g->Draw("AP");
    canvas->Draw("AP"); // Display the canvas object
    canvas->SaveAs(filename);

    delete canvas;
    delete g;
}
void create_graph(std::vector<double> x_values, std::vector<double> y_values,
                  char* title, char* x_label, char* y_label, char* filename) {
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
    g->GetXaxis()->SetTitle(x_label);
    g->GetYaxis()->SetTitle(y_label);
    g->SetMarkerStyle(kCircle);
    g->SetMarkerColor(kBlue);
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
    std::cout << "[";
    for (double i: v) {
        std::cout << i << ", ";
    }
    std::cout << "\n";
}

double get_var_m(double avg_m2, double avg_m) {
    double var = avg_m2 - std::pow(avg_m, 2);
    return var;
}

double onsager(double T) {
    double os = 0.0;
    if (T <= Tc) {
        os = std::pow(1 - std::pow(std::sinh(2/T), -4), 1/8);
    }
    return os;
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

double magnetic_susceptibility(double N, double T, double var_m) {
    //double var = avg_m2 - std::pow(avg_m, 2);
    double ms = N / T * var_m; // k considered k = 1
    return ms;
}

double magnetization_error(double var_m, double c) {
    double tau_m = 0.0;
    if (c != 1.0) {
        tau_m = c / (1.0 - c);
    }
    double m_err = std::sqrt(var_m * (2 * tau_m + 1) / NUM_MC);
    return m_err;
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

void ising_metropolis(int** s, int L, int N,  double T, double* avg_m, double* avg_m2, double* avg_m4, double* c,
                         std::mt19937 GEN,
                         std::uniform_real_distribution<> uniform, std::uniform_int_distribution<> select_point_uniform) {

    double m_aux = 1.0;
    for (int mc = 0; mc < NUM_MC; mc++) {
        // Montecarlo Step L*L
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

        double m =  std::abs(magnetization(s, L, N));
        *avg_m += m;
        *avg_m2 += std::pow(m, 2);
        *avg_m4 += std::pow(m, 4);

        *c = *c + m * m_aux;
        m_aux = m;
    }
    *avg_m /= NUM_MC;
    *avg_m2 /= NUM_MC;
    *avg_m4 /= NUM_MC;
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
        std::vector<double> temperatures;
        std::vector<double> onsager_values;
        std::vector<double> avg_m_values;
        std::vector<double> avg_m_err_values;
        std::vector<double> ms_values;
        std::vector<double> binder_values;

        for (int t = 0; t < NUM_T; t++) {
            const double T = TEMPERATURES[t];
            std::cout << "Computing T = " << T << std::endl;

            double avg_m = 0.0;
            double avg_m2 = 0.0;
            double avg_m4 = 0.0;
            double c = 0.0;

            unordered_initialize_s(s, L, GEN, uniform);
            ising_metropolis(s, L, N, T, &avg_m, &avg_m2, &avg_m4, &c, GEN, uniform, select_point_uniform);

            // Auxiliary calculations.
            double var_m = get_var_m(avg_m2, avg_m);

            // Compute the average magnetization error.
            c = (c / NUM_MC - std::pow(avg_m, 2)) / var_m;
            double avg_m_err = magnetization_error(var_m, c);

            // Compute the magnetic susceptibility value.
            double ms = magnetic_susceptibility(N, T, var_m);

            // Compute the binder ratio value.
            double br = binder_ratio(avg_m4, avg_m2);

            // Print status.
            std::cout << "T = " << T << ", avg_m(avg_m_err) = " << avg_m << "(" << avg_m_err << ")"
                        << ", ms = " << ms << ", br = " << br << std::endl;

            // Store results.
            temperatures.push_back(T);
            onsager_values.push_back(onsager(T));
            avg_m_values.push_back(avg_m);
            avg_m_err_values.push_back(avg_m_err);
            ms_values.push_back(ms);
            binder_values.push_back(br);
        }

        // Print results.
        std::cout << "L = " << L << ", N = " << N << ", NUM_MC = " << NUM_MC << ", NUM_T = " << NUM_T << std::endl;
        std::cout << "#Temperatures: " << std::endl;
        print_std_v(temperatures);
        std::cout << "#Onsager: " << std::endl;
        print_std_v(onsager_values);
        std::cout << "#Average Magnetization: " << std::endl;
        print_std_v(avg_m_values);
        std::cout << "#Average Magnetization Errors: " << std::endl;
        print_std_v(avg_m_err_values);
        std::cout << "#Magnetic Susceptibility: " << std::endl;
        print_std_v(ms_values);
        std::cout << "#Binder Ratio: " << std::endl;
        print_std_v(binder_values);

        // Create graphs with the results.
        create_graph_err(temperatures, avg_m_values, avg_m_err_values,
                         "Average Magnetization with errors", "T", "m_N", "avg_m_err.png");
        create_graph(temperatures, ms_values,
                     "Magnetic susceptibility", "T", "chi_N", "ms.png");
        create_graph(temperatures, binder_values,
                     "Binder Ratio", "T", "U^4", "br.png");

        // Deallocation of the 2D Lattice.
        for (int i = 0; i < L; i++) {
            delete[] s[i];
        }
        delete[] s;
    }
    return 0;
}

