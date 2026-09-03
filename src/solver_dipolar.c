#include "structures_nonspherical.h"
#include "facdes2Y.h"
#include "math_aux.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Forward declarations of closure functions
void closure_MSA_dipolar(double **c, double **eta, double *r, int n_points, double beta_mu2, double sigma);
void closure_LHNC_dipolar(double **c, double **h, double **eta, double *r, int n_points, double beta_mu2, double sigma);
void closure_QHNC_dipolar(double **c, double **h, double **eta, double *r, int n_points, double beta_mu2, double sigma);
void closure_RHNC_dipolar(double **c, double **h, double **eta, double *r, int n_points, double beta_mu2, double sigma, double *c_HS, double *h_HS);

/**
 * @brief Solves a small linear system A * x = b of dimension n (n <= 6)
 *        using Gaussian elimination with partial pivoting.
 */
static int solve_small_system(int n, double A[6][6], double b[6], double x[6]) {
    double initial_norm = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(A[i][j]) > initial_norm) initial_norm = fabs(A[i][j]);
        }
    }
    if (initial_norm < 1e-30) return 0;

    for (int i = 0; i < n; i++) {
        int max_row = i;
        double max_val = fabs(A[i][i]);
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > max_val) {
                max_val = fabs(A[k][i]);
                max_row = k;
            }
        }
        if (max_val < 1e-12 * initial_norm) return 0; // Scale-invariant relative check
        if (max_row != i) {
            for (int k = i; k < n; k++) {
                double tmp = A[i][k];
                A[i][k] = A[max_row][k];
                A[max_row][k] = tmp;
            }
            double tmp = b[i];
            b[i] = b[max_row];
            b[max_row] = tmp;
        }
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        if (fabs(A[i][i]) < 1e-12 * initial_norm) return 0;
        double sum = b[i];
        for (int j = i + 1; j < n; j++) {
            sum -= A[i][j] * x[j];
        }
        x[i] = sum / A[i][i];
    }
    return 1;
}

/**
 * @brief Computes the exact analytical Fourier/Hankel transform of the dipolar tail 1/r^3.
 *
 * Outside the hard core (r > sigma), c^{112}(r) = beta*mu^2 / r^3.
 * The exact integral \int_sigma^\infty (1/r) * j_2(k*r) dr evaluates analytically to j_1(k*sigma)/(k*sigma).
 * Therefore, C_{tail}^{112}(k) = -4*pi*beta*mu^2 * j_1(k*sigma) / (k*sigma).
 */
static inline double analytical_dipolar_tail_k(double k_val, double beta_mu2, double sigma) {
    double x = k_val * sigma;
    double j1_over_x;
    if (x < 1e-4) {
        j1_over_x = 1.0 / 3.0 - (x * x) / 30.0 + (x * x * x * x) / 840.0;
    } else {
        j1_over_x = (sin(x) - x * cos(x)) / (x * x * x);
    }
    return -4.0 * M_PI * beta_mu2 * j1_over_x;
}

/**
 * @brief Loads a 4-column analytical structure factor file (k, S000, S110, S112)
 * and linearly interpolates the curves onto the solver's discrete k-grid.
 */
static int load_analytical_sk(const char *filepath, const double *k_grid, int nodes,
                              double *S000_out, double *S110_out, double *S112_out) {
    FILE *f = fopen(filepath, "r");
    if (!f) {
        fprintf(stderr, "Error: Could not open init-sk file: %s\n", filepath);
        return -1;
    }

    int cap = 2048;
    double *k_file = malloc(cap * sizeof(double));
    double *s000_file = malloc(cap * sizeof(double));
    double *s110_file = malloc(cap * sizeof(double));
    double *s112_file = malloc(cap * sizeof(double));
    int count = 0;

    char line[512];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        double kv, s0, s1, s2;
        if (sscanf(line, "%lf %lf %lf %lf", &kv, &s0, &s1, &s2) == 4) {
            if (count >= cap) {
                cap *= 2;
                k_file = realloc(k_file, cap * sizeof(double));
                s000_file = realloc(s000_file, cap * sizeof(double));
                s110_file = realloc(s110_file, cap * sizeof(double));
                s112_file = realloc(s112_file, cap * sizeof(double));
            }
            k_file[count] = kv;
            s000_file[count] = s0;
            s110_file[count] = s1;
            s112_file[count] = s2;
            count++;
        }
    }
    fclose(f);

    if (count < 2) {
        fprintf(stderr, "Error: Insufficient data points in %s (found %d)\n", filepath, count);
        free(k_file); free(s000_file); free(s110_file); free(s112_file);
        return -1;
    }

    for (int i = 0; i < nodes; i++) {
        double target_k = k_grid[i];
        if (target_k <= k_file[0]) {
            S000_out[i] = s000_file[0];
            S110_out[i] = s110_file[0];
            S112_out[i] = s112_file[0];
        } else if (target_k >= k_file[count - 1]) {
            double k_max_f = k_file[count - 1];
            double decay = exp(-(target_k - k_max_f) / 10.0);
            S000_out[i] = 1.0 + (s000_file[count - 1] - 1.0) * decay;
            S110_out[i] = s110_file[count - 1] * decay;
            S112_out[i] = s112_file[count - 1] * decay;
        } else {
            int low = 0, high = count - 1;
            while (low <= high) {
                int mid = (low + high) / 2;
                if (k_file[mid] <= target_k) {
                    low = mid + 1;
                } else {
                    high = mid - 1;
                }
            }
            int idx = high;
            if (idx < 0) idx = 0;
            if (idx >= count - 1) idx = count - 2;
            double t = (target_k - k_file[idx]) / (k_file[idx + 1] - k_file[idx]);
            S000_out[i] = s000_file[idx] + t * (s000_file[idx + 1] - s000_file[idx]);
            S110_out[i] = s110_file[idx] + t * (s110_file[idx + 1] - s110_file[idx]);
            S112_out[i] = s112_file[idx] + t * (s112_file[idx + 1] - s112_file[idx]);
        }
    }

    free(k_file);
    free(s000_file);
    free(s110_file);
    free(s112_file);
    return 0;
}

/**
 * @brief Solves the OZ equation in k-space for Dipolar Hard Spheres.
 *
 * Uses the chi-basis decoupling (Blum/Wertheim) which reduces the coupled
 * 2x2 system {110, 112} to two independent scalar equations.
 *
 * From the Patey OZ equations:
 *   eta^000 = rho * C000 * H000    (000 decouples)
 *   eta^0   = (rho/3) * C^0 * H^0  (chi=0 mode, C^0 = C110 + 2*C112)
 *   eta^1   = -(rho/3) * C^1 * H^1 (chi=1 mode, C^1 = C110 - C112)
 *
 * Parity (-1)^chi yields:
 *   denom0 = 1 - (rho/3)*C0
 *   denom1 = 1 + (rho/3)*C1
 *
 * Recovery: H110 = (H^0 + 2*H^1)/3,  H112 = (H^0 - H^1)/3
 */
void solve_oz_k_space(double **C_k, double **H_k, int nodes, double rho) {
    for (int i = 0; i < nodes; i++) {
        double C000 = C_k[0][i];
        double C110 = C_k[1][i];
        double C112 = C_k[2][i];

        // --- 000 mode (decoupled) ---
        double denom000 = 1.0 - rho * C000;
        H_k[0][i] = (fabs(denom000) > 1e-12) ? C000 / denom000 : 0.0;

        // --- Chi modes ---
        double C0 = C110 + 2.0 * C112;  // chi=0: f^0 = f^{110} + 2*f^{112}
        double C1 = C110 - C112;         // chi=1: f^1 = f^{110} - f^{112}

        double denom0 = 1.0 - (rho / 3.0) * C0;
        double denom1 = 1.0 - (rho / 3.0) * C1; // Consistent with analytical MSA convention

        double H0 = (fabs(denom0) > 1e-12) ? C0 / denom0 : 0.0;
        double H1 = (fabs(denom1) > 1e-12) ? C1 / denom1 : 0.0;

        // --- Recover Patey projections ---
        H_k[1][i] = (H0 + 2.0 * H1) / 3.0;  // H110
        H_k[2][i] = (H0 - H1) / 3.0;         // H112
    }
}

/**
 * @brief Computes the exact Percus-Yevick Hard Sphere reference functions.
 */
void compute_HS_reference(double *c_HS, double *h_HS, double *r, double *k, 
                          int nodes, double dr, double rho, double sigma,
                          const double *K0_fwd, const double *K0_inv) {
    double eta_vol = rho * M_PI * pow(sigma, 3) / 6.0;
    double lambda1 = pow(1.0 + 2.0*eta_vol, 2) / pow(1.0 - eta_vol, 4);
    double lambda2 = -pow(1.0 + 0.5*eta_vol, 2) / pow(1.0 - eta_vol, 4);

    // 1. Exact PY c(r)
    for (int i = 0; i < nodes; i++) {
        if (r[i] < sigma) {
            c_HS[i] = -lambda1 - 6.0*eta_vol*lambda2*r[i] - 0.5*eta_vol*lambda1*pow(r[i], 3);
        } else {
            c_HS[i] = 0.0;
        }
    }

    // 2. Transform to C(k)
    double *C_k = malloc(nodes * sizeof(double));
    for (int i = 0; i < nodes; i++) {
        double sum = 0.0;
        const double *row = &K0_fwd[i * nodes];
        for (int j = 0; j < nodes; j++) {
            sum += row[j] * c_HS[j];
        }
        C_k[i] = sum;
    }

    // 3. Solve OZ for H(k)
    double *H_k = malloc(nodes * sizeof(double));
    for (int i = 0; i < nodes; i++) {
        double denom = 1.0 - rho * C_k[i];
        H_k[i] = (fabs(denom) > 1e-12) ? C_k[i] / denom : 0.0;
    }

    // 4. Transform back to h(r)
    for (int i = 0; i < nodes; i++) {
        double sum = 0.0;
        const double *row = &K0_inv[i * nodes];
        for (int j = 0; j < nodes; j++) {
            sum += row[j] * H_k[j];
        }
        h_HS[i] = sum;
    }

    free(C_k);
    free(H_k);
}

/**
 * @brief Main solver function for Dipolar Hard Spheres.
 */
void solver_dipolar(int closureID, double temp, double rho, double dipole_moment, 
                    int nodes, double rmax, const char *output_dir,
                    double temp_start, int ramp_steps, const char *init_sk_file) {
    
    printf("Initializing Dipolar Solver...\n");
    printf("Closure: %d (0=MSA, 1=LHNC, 2=QHNC, 3=RHNC)\n", closureID);
    printf("Dipole Moment (mu): %.4f\n", dipole_moment);
    printf("Density (rho): %.4f\n", rho);
    printf("rmax: %.2f, nodes: %d\n", rmax, nodes);

    // 1. Initialize Data Structures
    int n_projections = 3; // 0: 000, 1: 110, 2: 112
    ProjectionMatrix *h = create_projection_matrix(n_projections, nodes);
    ProjectionMatrix *c = create_projection_matrix(n_projections, nodes);
    ProjectionMatrix *eta = create_projection_matrix(n_projections, nodes);
    
    ProjectionMatrix *C_k = create_projection_matrix(n_projections, nodes);
    ProjectionMatrix *H_k = create_projection_matrix(n_projections, nodes);

    // Grid generation
    double dr = rmax / nodes;
    double *r = malloc(nodes * sizeof(double));
    double *k = malloc(nodes * sizeof(double));
    double dk = M_PI / rmax;

    for (int i = 0; i < nodes; i++) {
        r[i] = (i + 1) * dr;
        k[i] = (i + 1) * dk;
    }

    // Precompute Hankel transform kernel matrices for fast and orthogonal execution
    double *K0_fwd = malloc(nodes * nodes * sizeof(double));
    double *K0_inv = malloc(nodes * nodes * sizeof(double));
    double *K2_fwd = malloc(nodes * nodes * sizeof(double));
    double *K2_inv = malloc(nodes * nodes * sizeof(double));

    for (int i = 0; i < nodes; i++) {
        for (int j = 0; j < nodes; j++) {
            double kr = k[i] * r[j];
            double j0 = (fabs(kr) < 1e-7) ? 1.0 : sin(kr) / kr;
            double j2 = spherical_bessel_j2(kr);

            // Forward: C(k_i) = 4*pi*dr * sum_j r_j^2 * c(r_j) * j_l(k_i*r_j)
            K0_fwd[i * nodes + j] = 4.0 * M_PI * dr * (r[j] * r[j]) * j0;
            K2_fwd[i * nodes + j] = -4.0 * M_PI * dr * (r[j] * r[j]) * j2;

            // Inverse: h(r_i) = (dk / (2*pi^2)) * sum_j k_j^2 * H(k_j) * j_l(k_j*r_i)
            K0_inv[i * nodes + j] = (dk / (2.0 * M_PI * M_PI)) * (k[j] * k[j]) * j0;
            K2_inv[i * nodes + j] = -(dk / (2.0 * M_PI * M_PI)) * (k[j] * k[j]) * j2;
        }
    }

    double sigma = 1.0;

    // Hard Sphere Reference for RHNC
    double *c_HS = NULL;
    double *h_HS = NULL;
    if (closureID == 3) {
        c_HS = malloc(nodes * sizeof(double));
        h_HS = malloc(nodes * sizeof(double));
        compute_HS_reference(c_HS, h_HS, r, k, nodes, dr, rho, sigma, K0_fwd, K0_inv);
    }

    // Determine Temperature Continuation Schedule
    int num_stages = (init_sk_file != NULL && strlen(init_sk_file) > 0) ? 1 :
                     ((temp_start > temp && ramp_steps > 1) ? ramp_steps : 1);
    double *T_stages = malloc(num_stages * sizeof(double));
    if (num_stages > 1) {
        double ratio = pow(temp / temp_start, 1.0 / (num_stages - 1));
        for (int s = 0; s < num_stages; s++) {
            T_stages[s] = temp_start * pow(ratio, s);
        }
        T_stages[num_stages - 1] = temp;
        printf("Temperature Continuation Enabled: %d stages from T*=%.4f down to T*=%.4f\n", 
               num_stages, temp_start, temp);
    } else {
        T_stages[0] = temp;
    }

    // History for Anderson mixing (depth M = 4)
    const int M_DEPTH = 4;
    int total_dim = n_projections * nodes;
    double *c_hist[4];
    double *d_hist[4];
    for (int m = 0; m < M_DEPTH; m++) {
        c_hist[m] = malloc(total_dim * sizeof(double));
        d_hist[m] = malloc(total_dim * sizeof(double));
    }

    ProjectionMatrix *c_new_mat = create_projection_matrix(n_projections, nodes);
    double tolerance = 1e-6;

    // --- Temperature Continuation Loop ---
    for (int stage = 0; stage < num_stages; stage++) {
        double current_T = T_stages[stage];
        double beta = 1.0 / current_T;
        double beta_mu2 = beta * dipole_moment * dipole_moment;

        if (num_stages > 1) {
            printf("\n>>> [Continuation Stage %d/%d] Target T* = %.4f (beta*mu^2 = %.3f) <<<\n", 
                   stage + 1, num_stages, current_T, beta_mu2);
        }

        if (stage == 0) {
            if (init_sk_file != NULL && strlen(init_sk_file) > 0) {
                printf("Loading Analytical S(k) from %s for Direct Warm-Start...\n", init_sk_file);
                double *S000_in = malloc(nodes * sizeof(double));
                double *S110_in = malloc(nodes * sizeof(double));
                double *S112_in = malloc(nodes * sizeof(double));
                if (load_analytical_sk(init_sk_file, k, nodes, S000_in, S110_in, S112_in) == 0) {
                    double *C000_k = malloc(nodes * sizeof(double));
                    double *C110_k = malloc(nodes * sizeof(double));
                    double *C112_core_k = malloc(nodes * sizeof(double));

                    for (int i = 0; i < nodes; i++) {
                        C000_k[i] = (1.0 - 1.0 / S000_in[i]) / rho;

                        double S0_val = (S110_in[i] + 1.0) + 2.0 * S112_in[i];
                        double S1_val = (S110_in[i] + 1.0) - S112_in[i];

                        double C0_val = (3.0 / rho) * (1.0 - 1.0 / S0_val);
                        double C1_val = (3.0 / rho) * (1.0 - 1.0 / S1_val);

                        C110_k[i] = (C0_val + 2.0 * C1_val) / 3.0;
                        double C112_val = (C0_val - C1_val) / 3.0;

                        double C_tail_k_val = analytical_dipolar_tail_k(k[i], beta_mu2, sigma);
                        C112_core_k[i] = C112_val - C_tail_k_val;
                    }

                    for (int j = 0; j < nodes; j++) {
                        double c0_val = 0.0, c1_val = 0.0, c2_val = 0.0;
                        const double *row0 = &K0_inv[j * nodes];
                        const double *row2 = &K2_inv[j * nodes];
                        for (int i = 0; i < nodes; i++) {
                            c0_val += row0[i] * C000_k[i];
                            c1_val += row0[i] * C110_k[i];
                            c2_val += row2[i] * C112_core_k[i];
                        }
                        if (r[j] < sigma) {
                            c->data[0][j] = c0_val;
                            c->data[1][j] = c1_val;
                            c->data[2][j] = c2_val;
                        } else {
                            c->data[0][j] = 0.0;
                            c->data[1][j] = 0.0;
                            c->data[2][j] = beta_mu2 / pow(r[j], 3.0);
                        }
                    }

                    free(C000_k);
                    free(C110_k);
                    free(C112_core_k);
                    printf("Analytical direct correlation initialization c(r) complete.\n");
                } else {
                    closure_MSA_dipolar(c->data, eta->data, r, nodes, beta_mu2, sigma);
                }
                free(S000_in);
                free(S110_in);
                free(S112_in);
            } else {
                closure_MSA_dipolar(c->data, eta->data, r, nodes, beta_mu2, sigma);
            }
        } else {
            // Update dipole tail outside core for the new temperature
            for (int i = 0; i < nodes; i++) {
                if (r[i] > sigma) {
                    c->data[2][i] = beta_mu2 / pow(r[i], 3.0);
                }
            }
        }

        double alpha = (current_T < 0.3) ? 0.2 : 0.4; // Adaptive Anderson/Picard damping factor
        int hist_count = 0;
        int iter = 0;
        double error = 1.0;
        double stage_tol = (stage < num_stages - 1) ? 1e-4 : tolerance;
        int max_iter_stage = (stage < num_stages - 1) ? 300 : 600;

        while (iter < max_iter_stage && error > stage_tol) {
            
            // A. Forward Hankel Transforms c(r) -> C(k) with Dipolar Tail Splitting
            for (int i = 0; i < nodes; i++) {
                double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
                const double *row0 = &K0_fwd[i * nodes];
                const double *row2 = &K2_fwd[i * nodes];
                for (int j = 0; j < nodes; j++) {
                    sum0 += row0[j] * c->data[0][j];
                    sum1 += row0[j] * c->data[1][j];
                    // Core direct correlation function (strictly zero outside sigma in MSA)
                    double c_core_2 = (r[j] < sigma) ? c->data[2][j] : (c->data[2][j] - beta_mu2 / pow(r[j], 3.0));
                    sum2 += row2[j] * c_core_2;
                }
                C_k->data[0][i] = sum0;
                C_k->data[1][i] = sum1;
                // Add exact analytical Hankel transform of long-range 1/r^3 tail
                C_k->data[2][i] = sum2 + analytical_dipolar_tail_k(k[i], beta_mu2, sigma);
            }

            // B. Solve OZ in k-space
            solve_oz_k_space(C_k->data, H_k->data, nodes, rho);

            // C. Inverse Hankel Transforms N(k) = H(k) - C(k) -> eta(r)
            for (int i = 0; i < nodes; i++) {
                double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
                const double *row0 = &K0_inv[i * nodes];
                const double *row2 = &K2_inv[i * nodes];
                for (int j = 0; j < nodes; j++) {
                    double N0 = H_k->data[0][j] - C_k->data[0][j];
                    double N1 = H_k->data[1][j] - C_k->data[1][j];
                    double N2 = H_k->data[2][j] - C_k->data[2][j];
                    sum0 += row0[j] * N0;
                    sum1 += row0[j] * N1;
                    sum2 += row2[j] * N2;
                }
                eta->data[0][i] = sum0;
                eta->data[1][i] = sum1;
                eta->data[2][i] = sum2;

                h->data[0][i] = c->data[0][i] + sum0;
                h->data[1][i] = c->data[1][i] + sum1;
                h->data[2][i] = c->data[2][i] + sum2;
            }

            // D. Prepare target c_new buffer
            for (int p = 0; p < n_projections; p++) {
                for (int i = 0; i < nodes; i++) {
                    c_new_mat->data[p][i] = c->data[p][i];
                }
            }

            // E. Apply Closure Relation
            if (closureID == 0) {
                closure_MSA_dipolar(c_new_mat->data, eta->data, r, nodes, beta_mu2, sigma);
            } else if (closureID == 1) {
                closure_LHNC_dipolar(c_new_mat->data, h->data, eta->data, r, nodes, beta_mu2, sigma);
            } else if (closureID == 2) {
                closure_QHNC_dipolar(c_new_mat->data, h->data, eta->data, r, nodes, beta_mu2, sigma);
            } else if (closureID == 3) {
                closure_RHNC_dipolar(c_new_mat->data, h->data, eta->data, r, nodes, beta_mu2, sigma, c_HS, h_HS);
            }

            // F. Compute L2 residual error d = c_new - c
            error = 0.0;
            double *d_curr = malloc(total_dim * sizeof(double));
            double *c_curr = malloc(total_dim * sizeof(double));
            int idx = 0;
            for (int p = 0; p < n_projections; p++) {
                for (int i = 0; i < nodes; i++) {
                    double diff = c_new_mat->data[p][i] - c->data[p][i];
                    d_curr[idx] = diff;
                    c_curr[idx] = c->data[p][i];
                    error += diff * diff;
                    idx++;
                }
            }
            error = sqrt(error / total_dim);

            if (iter % 10 == 0 || error <= stage_tol) {
                printf("Iter %4d: Error = %.5e\n", iter, error);
            }

            if (error <= stage_tol) {
                free(d_curr);
                free(c_curr);
                iter++;
                break;
            }

            // Update Anderson history buffer
            if (hist_count < M_DEPTH) {
                memcpy(c_hist[hist_count], c_curr, total_dim * sizeof(double));
                memcpy(d_hist[hist_count], d_curr, total_dim * sizeof(double));
                hist_count++;
            } else {
                // Shift history
                double *tmp_c = c_hist[0];
                double *tmp_d = d_hist[0];
                for (int m = 0; m < M_DEPTH - 1; m++) {
                    c_hist[m] = c_hist[m + 1];
                    d_hist[m] = d_hist[m + 1];
                }
                c_hist[M_DEPTH - 1] = tmp_c;
                d_hist[M_DEPTH - 1] = tmp_d;
                memcpy(c_hist[M_DEPTH - 1], c_curr, total_dim * sizeof(double));
                memcpy(d_hist[M_DEPTH - 1], d_curr, total_dim * sizeof(double));
            }

            // Apply Anderson Mixing if history >= 2 (Walker-Ni formulation)
            int applied_anderson = 0;
            if (hist_count >= 2) {
                int m = hist_count;
                int K = m - 1; // Dimension of difference system
                double Amat[6][6] = {0};
                double bvec[6] = {0};
                double gamma_vec[6] = {0};

                // Compute Delta d_j = d_hist[j] - d_hist[K]
                // Amat[j][k] = <Delta d_j, Delta d_k>, bvec[j] = -<Delta d_j, d_hist[K]>
                for (int j = 0; j < K; j++) {
                    for (int k_idx = j; k_idx < K; k_idx++) {
                        double dot = 0.0;
                        for (int idx_pt = 0; idx_pt < total_dim; idx_pt++) {
                            double dd_j = d_hist[j][idx_pt] - d_hist[K][idx_pt];
                            double dd_k = d_hist[k_idx][idx_pt] - d_hist[K][idx_pt];
                            dot += dd_j * dd_k;
                        }
                        Amat[j][k_idx] = dot;
                        Amat[k_idx][j] = dot;
                    }

                    double dot_rhs = 0.0;
                    for (int idx_pt = 0; idx_pt < total_dim; idx_pt++) {
                        double dd_j = d_hist[j][idx_pt] - d_hist[K][idx_pt];
                        dot_rhs += dd_j * d_hist[K][idx_pt];
                    }
                    bvec[j] = -dot_rhs;
                }

                // Scale-invariant Tikhonov regularization
                double max_diag = 0.0;
                for (int j = 0; j < K; j++) {
                    if (Amat[j][j] > max_diag) max_diag = Amat[j][j];
                }
                double lambda = (max_diag > 1e-16) ? 1e-4 * max_diag : 1e-12;
                for (int j = 0; j < K; j++) {
                    Amat[j][j] += lambda;
                }

                if (solve_small_system(K, Amat, bvec, gamma_vec)) {
                    double gamma_last = 1.0;
                    for (int j = 0; j < K; j++) {
                        gamma_last -= gamma_vec[j];
                    }

                    // Smoothly bound extrapolation weights if excessive
                    double max_w = fabs(gamma_last);
                    for (int j = 0; j < K; j++) {
                        if (fabs(gamma_vec[j]) > max_w) max_w = fabs(gamma_vec[j]);
                    }
                    if (max_w > 5.0) {
                        double scale = 5.0 / max_w;
                        for (int j = 0; j < K; j++) gamma_vec[j] *= scale;
                        gamma_last = 1.0;
                        for (int j = 0; j < K; j++) gamma_last -= gamma_vec[j];
                    }

                    int val_valid = 1;
                    for (int pt = 0; pt < total_dim; pt++) {
                        double val = 0.0;
                        for (int j = 0; j < K; j++) {
                            val += gamma_vec[j] * (c_hist[j][pt] + alpha * d_hist[j][pt]);
                        }
                        val += gamma_last * (c_hist[K][pt] + alpha * d_hist[K][pt]);
                        if (isnan(val) || isinf(val) || fabs(val) > 1e6) {
                            val_valid = 0;
                            break;
                        }
                    }

                    if (val_valid) {
                        idx = 0;
                        for (int p = 0; p < n_projections; p++) {
                            for (int i = 0; i < nodes; i++) {
                                double val = 0.0;
                                for (int j = 0; j < K; j++) {
                                    val += gamma_vec[j] * (c_hist[j][idx] + alpha * d_hist[j][idx]);
                                }
                                val += gamma_last * (c_hist[K][idx] + alpha * d_hist[K][idx]);
                                c->data[p][i] = val;
                                idx++;
                            }
                        }
                        applied_anderson = 1;
                    } else {
                        hist_count = 0; // Discard invalid history
                    }
                }
            }

            // Fallback to standard damped Picard update
            if (!applied_anderson) {
                for (int p = 0; p < n_projections; p++) {
                    for (int i = 0; i < nodes; i++) {
                        c->data[p][i] += alpha * (c_new_mat->data[p][i] - c->data[p][i]);
                    }
                }
            }

            // In MSA closure, strictly enforce exact asymptotic behavior outside core
            if (closureID == 0) {
                for (int i = 0; i < nodes; i++) {
                    if (r[i] >= sigma) {
                        c->data[0][i] = 0.0;
                        c->data[1][i] = 0.0;
                        c->data[2][i] = beta_mu2 / pow(r[i], 3.0);
                    }
                }
            }

            free(d_curr);
            free(c_curr);
            iter++;
        }

        if (num_stages > 1) {
            printf(">>> Stage %d/%d (T*=%.4f) completed in %d iterations (Error = %.3e) <<<\n", 
                   stage + 1, num_stages, current_T, iter, error);
        }
    }

    free(T_stages);
    free_projection_matrix(c_new_mat);

    // Free Anderson history
    for (int m = 0; m < M_DEPTH; m++) {
        free(c_hist[m]);
        free(d_hist[m]);
    }

    // Final exact calculation of C_k and H_k from converged c(r)
    double final_beta_mu2 = (dipole_moment * dipole_moment) / temp;
    for (int i = 0; i < nodes; i++) {
        double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
        const double *row0 = &K0_fwd[i * nodes];
        const double *row2 = &K2_fwd[i * nodes];
        for (int j = 0; j < nodes; j++) {
            sum0 += row0[j] * c->data[0][j];
            sum1 += row0[j] * c->data[1][j];
            double c_core_2 = (r[j] < sigma) ? c->data[2][j] : (c->data[2][j] - final_beta_mu2 / pow(r[j], 3.0));
            sum2 += row2[j] * c_core_2;
        }
        C_k->data[0][i] = sum0;
        C_k->data[1][i] = sum1;
        C_k->data[2][i] = sum2 + analytical_dipolar_tail_k(k[i], final_beta_mu2, sigma);
    }
    solve_oz_k_space(C_k->data, H_k->data, nodes, rho);

    // 4. Output Results
    FILE *fp = fopen("output/output_dipolar.dat", "w");
    if (fp) {
        fprintf(fp, "# r h000 h110 h112 c000 c110 c112\n");
        for (int i = 0; i < nodes; i++) {
            fprintf(fp, "%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", 
                r[i], 
                h->data[0][i], h->data[1][i], h->data[2][i],
                c->data[0][i], c->data[1][i], c->data[2][i]);
        }
        fclose(fp);
        printf("Written output/output_dipolar.dat\n");
    }

    // Output k-space Results
    FILE *fp_k = fopen("output/output_dipolar_k.dat", "w");
    if (fp_k) {
        fprintf(fp_k, "# k H000 H110 H112 C000 C110 C112\n");
        for (int i = 0; i < nodes; i++) {
            fprintf(fp_k, "%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", 
                k[i], 
                H_k->data[0][i], H_k->data[1][i], H_k->data[2][i],
                C_k->data[0][i], C_k->data[1][i], C_k->data[2][i]);
        }
        fclose(fp_k);
        printf("Written output/output_dipolar_k.dat\n");
    }

    // ---------------------------------------------------------------
    // Output: Structure Factor S(k) in Patey and chi representations
    // ---------------------------------------------------------------
    FILE *fp_sk = fopen("output/output_dipolar_sk.dat", "w");
    if (fp_sk) {
        fprintf(fp_sk, "# k S000 S110 S112 S10 S11\n");
        for (int i = 0; i < nodes; i++) {
            double C110 = C_k->data[1][i];
            double C112 = C_k->data[2][i];

            // --- Chi representation (from C_k) ---
            double C0 = C110 + 2.0 * C112;  // chi=0 mode
            double C1 = C110 - C112;         // chi=1 mode

            double denom0 = 1.0 - (rho / 3.0) * C0;
            double denom1 = 1.0 - (rho / 3.0) * C1;

            double S0 = (fabs(denom0) > 1e-12) ? 1.0 / denom0 : 1e12;
            double S1 = (fabs(denom1) > 1e-12) ? 1.0 / denom1 : 1e12;

            double S000 = 1.0 + rho * H_k->data[0][i];
            double S110 = (S0 + 2.0 * S1) / 3.0 - 1.0;
            double S112 = (S0 - S1) / 3.0;

            fprintf(fp_sk, "%.5e  %.5e  %.5e  %.5e  %.5e  %.5e\n",
                k[i], S000, S110, S112, S0, S1);
        }
        fclose(fp_sk);
        printf("Written output/output_dipolar_sk.dat\n");
    }

    // Cleanup
    free_projection_matrix(h);
    free_projection_matrix(c);
    free_projection_matrix(eta);
    free_projection_matrix(C_k);
    free_projection_matrix(H_k);
    free(r);
    free(k);
    free(K0_fwd);
    free(K0_inv);
    free(K2_fwd);
    free(K2_inv);
    if (c_HS) free(c_HS);
    if (h_HS) free(h_HS);
}
