#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

struct Options {
    std::string density = "p1";
    int64_t samples = 20000;
    int64_t burn_in = 5000;
    int64_t thin = 1;
    double step = 0.6;
    uint64_t seed = 12345;
    int grid_nx = 200;
    int grid_ny = 200;
    std::string output = "";
};

bool in_domain(double x, double y) {
    if (x < -2.0 || x > 2.0) {
        return false;
    }
    if (y < 0.0) {
        return false;
    }
    double ymax = 2.0 + x;
    return y <= ymax;
}

double log_q1(double x, double y) {
    return -(x * x) - (y * y) + x * y;
}

double log_q2(double x, double y) {
    return -(x * x) - (y * y);
}

double simpson_1d(const std::vector<double> &values, double h) {
    if (values.size() < 2) {
        return 0.0;
    }
    double sum = values.front() + values.back();
    for (size_t i = 1; i + 1 < values.size(); ++i) {
        sum += (i % 2 == 0) ? 2.0 * values[i] : 4.0 * values[i];
    }
    return sum * h / 3.0;
}

double integrate_2d_simpson(int nx, int ny, const std::string &density) {
    if (nx % 2 == 1) {
        nx += 1;
    }
    if (ny % 2 == 1) {
        ny += 1;
    }
    double x_min = -2.0;
    double x_max = 2.0;
    double hx = (x_max - x_min) / nx;
    std::vector<double> x_values(nx + 1);
    for (int i = 0; i <= nx; ++i) {
        x_values[i] = x_min + i * hx;
    }

    std::vector<double> fx(nx + 1, 0.0);
    for (int i = 0; i <= nx; ++i) {
        double x = x_values[i];
        double y_min = 0.0;
        double y_max = 2.0 + x;
        if (y_max <= y_min) {
            fx[i] = 0.0;
            continue;
        }
        double hy = (y_max - y_min) / ny;
        std::vector<double> y_vals(ny + 1, 0.0);
        for (int j = 0; j <= ny; ++j) {
            double y = y_min + j * hy;
            double val = 0.0;
            if (density == "p1") {
                val = std::exp(log_q1(x, y));
            } else {
                val = std::exp(log_q2(x, y));
            }
            y_vals[j] = val;
        }
        fx[i] = simpson_1d(y_vals, hy);
    }
    return simpson_1d(fx, hx);
}

Options parse_args(int argc, char **argv) {
    Options opts;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto next = [&](double &target) {
            if (i + 1 < argc) {
                target = std::stod(argv[++i]);
            }
        };
        auto next_int = [&](int64_t &target) {
            if (i + 1 < argc) {
                target = std::stoll(argv[++i]);
            }
        };
        auto next_uint = [&](uint64_t &target) {
            if (i + 1 < argc) {
                target = static_cast<uint64_t>(std::stoull(argv[++i]));
            }
        };
        auto next_str = [&](std::string &target) {
            if (i + 1 < argc) {
                target = argv[++i];
            }
        };

        if (arg == "--density") {
            next_str(opts.density);
        } else if (arg == "--samples") {
            next_int(opts.samples);
        } else if (arg == "--burn") {
            next_int(opts.burn_in);
        } else if (arg == "--thin") {
            next_int(opts.thin);
        } else if (arg == "--step") {
            next(opts.step);
        } else if (arg == "--seed") {
            next_uint(opts.seed);
        } else if (arg == "--grid-nx") {
            int64_t tmp = opts.grid_nx;
            next_int(tmp);
            opts.grid_nx = static_cast<int>(tmp);
        } else if (arg == "--grid-ny") {
            int64_t tmp = opts.grid_ny;
            next_int(tmp);
            opts.grid_ny = static_cast<int>(tmp);
        } else if (arg == "--output") {
            next_str(opts.output);
        }
    }
    return opts;
}

int main(int argc, char **argv) {
    Options opts = parse_args(argc, argv);
    if (opts.density != "p1" && opts.density != "p2") {
        std::cerr << "Unknown density: " << opts.density << "\n";
        return 1;
    }

    if (opts.samples <= 0 || opts.burn_in < 0 || opts.thin <= 0) {
        std::cerr << "Invalid sampling parameters.\n";
        return 1;
    }

    double z_norm = integrate_2d_simpson(opts.grid_nx, opts.grid_ny, opts.density);

    std::mt19937_64 rng(opts.seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    double x = 0.0;
    double y = 1.0;
    if (!in_domain(x, y)) {
        x = 0.0;
        y = 0.0;
    }

    auto log_q = [&](double x_val, double y_val) {
        return (opts.density == "p1") ? log_q1(x_val, y_val) : log_q2(x_val, y_val);
    };

    int64_t total_steps = opts.burn_in + opts.samples * opts.thin;
    int64_t accepted = 0;
    int64_t kept = 0;

    std::ofstream out;
    if (!opts.output.empty()) {
        out.open(opts.output);
        out << "idx,x,y,observable\n";
    }

    double mean = 0.0;
    double m2 = 0.0;

    for (int64_t step_idx = 0; step_idx < total_steps; ++step_idx) {
        double x_prop = x + opts.step * normal(rng);
        double y_prop = y + opts.step * normal(rng);
        bool accept = false;
        if (in_domain(x_prop, y_prop)) {
            double log_alpha = log_q(x_prop, y_prop) - log_q(x, y);
            if (std::log(uniform(rng)) < log_alpha) {
                accept = true;
            }
        }
        if (accept) {
            x = x_prop;
            y = y_prop;
            accepted++;
        }

        if (step_idx >= opts.burn_in && ((step_idx - opts.burn_in) % opts.thin == 0)) {
            double h = 2.0 + x * x + y * y;
            double observable = h;
            if (opts.density == "p2") {
                observable = h * std::exp(x * y);
            }
            kept++;
            double delta = observable - mean;
            mean += delta / kept;
            m2 += delta * (observable - mean);

            if (out.is_open()) {
                out << kept << "," << std::setprecision(10) << x << "," << y << "," << observable << "\n";
            }
        }
    }

    double variance = (kept > 1) ? (m2 / (kept - 1)) : 0.0;
    double stddev = std::sqrt(variance);
    double estimate = z_norm * mean;
    double stderr = (kept > 0) ? z_norm * stddev / std::sqrt(static_cast<double>(kept)) : 0.0;

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Density: " << opts.density << "\n";
    std::cout << "Z_norm: " << z_norm << "\n";
    std::cout << "Samples kept: " << kept << "\n";
    std::cout << "Acceptance rate: " << static_cast<double>(accepted) / total_steps << "\n";
    std::cout << "Mean observable: " << mean << "\n";
    std::cout << "Integral estimate: " << estimate << " +/- " << stderr << "\n";

    return 0;
}
