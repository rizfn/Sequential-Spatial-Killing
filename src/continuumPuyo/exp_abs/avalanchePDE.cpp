#include <random>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cmath>
#include <thread>
#include <mutex>
#include <algorithm>

#pragma GCC optimize("O3", "inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);
    return 0;
}();

// High stability, fast constants
constexpr int L = 128;
constexpr double T = 8000.0;
constexpr double dt = 0.01;
constexpr double a = 0.1;
constexpr double c_val = 1;
constexpr double lambda_exp = 100;
constexpr double lambda_abs = 1;
constexpr double dx = 0.1;
constexpr double noise_strength = 0.1;
constexpr int num_runs = 12;
constexpr int record_steps = 10000;
constexpr double hist_min = -2.5;
constexpr double hist_max = 2.5;
constexpr int hist_bins = 400;

void run_single_sim(int seed, int N, std::vector<double>& out_roughness, std::vector<double>& out_profile, std::vector<uint64_t>& out_slope_hist) {
    // Fast linear congruential generator (much faster than mt19937 for massive loops)
    std::minstd_rand gen(seed);
    std::normal_distribution<double> dist(0.0, 1.0);

    long long steps = std::round(T / dt);
    
    double p = 1.0 / N;
    double c_diff = (a * a * p * p) / 2.0;
    double c_nonlin_exp = 2.0 * lambda_exp / N; 
    double c_nonlin_abs = lambda_abs * a * p * p;
    
    // Spatiotemporal white noise must scale with 1/sqrt(dx)
    double sqrt_dt_dx = std::sqrt(dt / dx);

    // Number of array elements is Macroscopic Size / grid scaling
    int array_size = std::round(L / dx);

    std::vector<double> h(array_size, 0.0);
    std::vector<double> new_h(array_size, 0.0);
    std::vector<double> k1(array_size, 0.0);
    std::vector<double> k2(array_size, 0.0);
    std::vector<double> h_tmp(array_size, 0.0);
    std::vector<double> noise(array_size, 0.0);

    long long tail_steps = std::max(1LL, (long long)(steps * 0.1));
    double hist_bin_width = (hist_max - hist_min) / hist_bins;

    // We still want 'record_steps' total data points over the macroscopic time T
    long long record_interval = std::max(1LL, steps / record_steps);
    int record_idx = 0;

    // Precalculate division constants for massive speedup in the hot loop
    double inv_dx = 1.0 / dx;
    double inv_2dx = 1.0 / (2.0 * dx);
    double inv_3 = 1.0 / 3.0;

    auto calc_explicit = [&](const std::vector<double>& current_h, std::vector<double>& out_explicit) {
        for (int x = 0; x < array_size; ++x) {
            int left = (x == 0) ? array_size - 1 : x - 1;
            int right = (x == array_size - 1) ? 0 : x + 1;

            double grad_h_central = (current_h[right] - current_h[left]) * inv_2dx;
            double abs_grad_h = std::abs(grad_h_central);

            out_explicit[x] = c_nonlin_exp * std::exp(-c_val * abs_grad_h) - c_nonlin_abs * abs_grad_h;
        }
    };

    // Thomas algorithm for tridiagonal system augmented with Sherman-Morrison for cyclic boundaries
    auto solve_cyclic_tridiagonal = [&](double a, double b, double c, const std::vector<double>& rhs, std::vector<double>& out) {
        int n = array_size;
        std::vector<double> cp(n, 0.0);
        std::vector<double> dp1(n, 0.0);
        std::vector<double> dp2(n, 0.0);
        double gamma = -b;

        // Modified system 1: standard tridiagonal with b' = b - gamma at 0, b' = b - a*c/gamma at n-1
        double bb = b - gamma;
        cp[0] = c / bb;
        dp1[0] = rhs[0] / bb;
        dp2[0] = gamma / bb;

        for (int i = 1; i < n - 1; i++) {
            double m = 1.0 / (b - a * cp[i - 1]);
            cp[i] = c * m;
            dp1[i] = (rhs[i] - a * dp1[i - 1]) * m;
            dp2[i] = (0.0 - a * dp2[i - 1]) * m;
        }

        double b_last = b - a * c / gamma;
        double m = 1.0 / (b_last - a * cp[n - 2]);
        dp1[n - 1] = (rhs[n - 1] - a * dp1[n - 2]) * m;
        dp2[n - 1] = (c - a * dp2[n - 2]) * m;

        std::vector<double> y(n, 0.0);
        std::vector<double> q(n, 0.0);
        y[n - 1] = dp1[n - 1];
        q[n - 1] = dp2[n - 1];

        for (int i = n - 2; i >= 0; i--) {
            y[i] = dp1[i] - cp[i] * y[i + 1];
            q[i] = dp2[i] - cp[i] * q[i + 1];
        }

        double num = 0.0;
        double den = 1.0;
        num = y[0] + a / gamma * y[n - 1];
        den = 1.0 + q[0] + a / gamma * q[n - 1];

        double v_dot_y_over_1_plus_v_dot_q = num / den;

        for (int i = 0; i < n; i++) {
            out[i] = y[i] - q[i] * v_dot_y_over_1_plus_v_dot_q;
        }
    };

    double r = c_diff * dt * inv_dx * inv_dx;
    double diag_b = 1.0 + 2.0 * r;
    double offdiag_a = -r;
    double offdiag_c = -r;

    std::vector<double> h_star(array_size, 0.0);

    // Crank-Nicolson-like IMEX RK-SSP2 for stability of non-linear terms
    for (long long i = 0; i < steps; ++i) {
        // Stage 1: Explicit Euler predictor + Implicit solve
        calc_explicit(h, k1);
        for (int x = 0; x < array_size; ++x) {
            h_tmp[x] = h[x] + k1[x] * dt;
        }
        solve_cyclic_tridiagonal(offdiag_a, diag_b, offdiag_c, h_tmp, h_star);

        // Stage 2: 2nd order corrector + Implicit solve
        calc_explicit(h_star, k2);
        for (int x = 0; x < array_size; ++x) {
            noise[x] = noise_strength * dist(gen) * sqrt_dt_dx;
            // The un-inverted equation: 0.5 * h_old + 0.5 * (h_star + dt * k2) + noise
            h_tmp[x] = 0.5 * h[x] + 0.5 * h_star[x] + 0.5 * k2[x] * dt + noise[x];
        }
        
        // The implicit operator here needs to apply to the 0.5 * h_star component exactly like the system requires
        // Since the operator is linear, we can just solve the standard step using half off-diagonals, 
        // but for an SSP2 IMEX, applying the same implicit operator solve on the averaged explicit part is standard:
        solve_cyclic_tridiagonal(offdiag_a, diag_b, offdiag_c, h_tmp, new_h);

        for (int x = 0; x < array_size; ++x) {
            int left = (x == 0) ? array_size - 1 : x - 1;
            int right = (x == array_size - 1) ? 0 : x + 1;

            // Histogram local slopes at the tail
            if (i >= steps - tail_steps) {
                // For observables, central difference is still mathematically unbiased
                double grad_h_central = (new_h[right] - new_h[left]) / (2.0 * dx);
                int bin = std::floor((grad_h_central - hist_min) / hist_bin_width);
                if (bin >= 0 && bin < hist_bins) {
                    out_slope_hist[bin]++;
                }
            }
        }

        // Subsample roughness only when needed
        if (i % record_interval == 0 && record_idx < record_steps) {
            double mean_h = 0.0;
            for(int j=0; j<array_size; ++j) mean_h += h[j];
            mean_h /= array_size;

            double sq_diff = 0.0;
            for (int x = 0; x < array_size; ++x) {
                double diff = h[x] - mean_h;
                sq_diff += diff * diff;
            }

            out_roughness[record_idx] = std::sqrt(sq_diff / array_size);
            record_idx++;
        }

        std::swap(h, new_h);
        
        // Prevent blowups from destroying entire ensemble mid-run
        if (std::isnan(h[0])) {
            break; 
        }
    }

    double final_mean = 0.0;
    for(int j=0; j<array_size; ++j) final_mean += h[j];
    final_mean /= array_size;
    for(int j=0; j<array_size; ++j) out_profile[j] = h[j] - final_mean;
}

void simulate_ensemble(int N, const std::string& out_dir) {
    int array_size = std::round(L / dx);
    std::vector<std::vector<double>> all_roughness(num_runs, std::vector<double>(record_steps, 0.0));
    std::vector<std::vector<double>> all_profiles(num_runs, std::vector<double>(array_size, 0.0));
    std::vector<std::vector<uint64_t>> all_slope_hists(num_runs, std::vector<uint64_t>(hist_bins, 0));

    std::vector<std::thread> threads;
    for(int r=0; r<num_runs; ++r) {
        int seed = 42000 + r * 100 + N;
        threads.emplace_back(run_single_sim, seed, N, std::ref(all_roughness[r]), std::ref(all_profiles[r]), std::ref(all_slope_hists[r]));
    }

    for(auto& t : threads) {
        if(t.joinable()) t.join();
    }

    std::vector<double> avg_roughness(record_steps, 0.0);
    for(int r=0; r<num_runs; ++r) {
        for(int i=0; i<record_steps; ++i) {
            avg_roughness[i] += all_roughness[r][i] / num_runs;
        }
    }

    std::vector<uint64_t> combined_hist(hist_bins, 0);
    for(int r=0; r<num_runs; ++r) {
        for(int b=0; b<hist_bins; ++b) {
            combined_hist[b] += all_slope_hists[r][b];
        }
    }

    std::filesystem::path dir(out_dir);
    std::filesystem::create_directories(dir);

    // Save Roughness
    std::ofstream f_rough(dir / ("roughness_N" + std::to_string(N) + ".tsv"));
    f_rough << "time\troughness\n";
    for(int i=0; i<record_steps; ++i) {
        double time_val = i * (T / record_steps);
        f_rough << std::fixed << std::setprecision(6) << time_val << "\t" << avg_roughness[i] << "\n";
    }
    f_rough.close();

    // Save Profile (just take run 0 as rep)
    std::ofstream f_prof(dir / ("profile_N" + std::to_string(N) + ".tsv"));
    f_prof << "x\th\n";
    for(int x=0; x<array_size; ++x) {
        f_prof << std::fixed << std::setprecision(6) << x * dx << "\t" << all_profiles[0][x] << "\n";
    }
    f_prof.close();

    // Save Slopes
    std::ofstream f_slope(dir / ("slopes_N" + std::to_string(N) + ".tsv"));
    f_slope << "slope_bin_center\tcount\n";
    double bin_width = (hist_max - hist_min) / hist_bins;
    for(int b=0; b<hist_bins; ++b) {
        double center = hist_min + (b + 0.5) * bin_width;
        f_slope << std::fixed << std::setprecision(6) << center << "\t" << combined_hist[b] << "\n";
    }
    f_slope.close();
}

int main() {
    std::vector<int> Ns = {1, 2, 4, 6, 8, 12, 15, 20, 40};

    std::ostringstream oss;
    oss << "outputs/avalanchePDE_cpp/L" << L << "_T" << (int)T 
        << "_dt" << dt << "_a" << a << "_c" << c_val << "_lamexp" << lambda_exp << "_lam" << lambda_abs;
    std::string out_dir = oss.str();

    for(int N : Ns) {
        std::cout << "Running C++ ensemble for N=" << N << "..." << std::endl;
        simulate_ensemble(N, out_dir);
    }

    std::cout << "All simulations complete. Outputs in " << out_dir << std::endl;
    return 0;
}