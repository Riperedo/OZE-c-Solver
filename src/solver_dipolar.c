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
    for (int i = 0; i < n; i++) {
        int max_row = i;
        double max_val = fabs(A[i][i]);
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > max_val) {
                max_val = fabs(A[k][i]);
                max_row = k;
            }
        }
        if (max_val < 1e-13) return 0; // Singular matrix
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
        if (fabs(A[i][i]) < 1e-13) return 0;
        double sum = b[i];
        for (int j = i + 1; j < n; j++) {
            sum -= A[i][j] * x[j];
        }
        x[i] = sum / A[i][i];
    }
    return 1;
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
                   int nodes, double rmax, const char *output_dir) {
    
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

    double beta = 1.0 / temp;
    double beta_mu2 = beta * dipole_moment * dipole_moment;
    double sigma = 1.0;

    // Hard Sphere Reference for RHNC
    double *c_HS = NULL;
    double *h_HS = NULL;
    if (closureID == 3) {
        c_HS = malloc(nodes * sizeof(double));
        h_HS = malloc(nodes * sizeof(double));
        compute_HS_reference(c_HS, h_HS, r, k, nodes, dr, rho, sigma, K0_fwd, K0_inv);
    }

    // 2. Initialization 
    closure_MSA_dipolar(c->data, eta->data, r, nodes, beta_mu2, sigma);

    // 3. Iteration Loop with Anderson Mixing Acceleration
    int max_iter = 1000;
    double tolerance = 1e-6;
    double error = 1.0;
    int iter = 0;
    double alpha = 0.5; // Anderson damping factor

    // History for Anderson mixing (depth M = 4)
    const int M_DEPTH = 4;
    int hist_count = 0;
    int total_dim = n_projections * nodes;
    double *c_hist[4];
    double *d_hist[4];
    for (int m = 0; m < M_DEPTH; m++) {
        c_hist[m] = malloc(total_dim * sizeof(double));
        d_hist[m] = malloc(total_dim * sizeof(double));
    }

    printf("Starting Iteration with Anderson Mixing Acceleration...\n");
    ProjectionMatrix *c_new_mat = create_projection_matrix(n_projections, nodes);

    while (iter < max_iter && error > tolerance) {
        
        // A. Forward Hankel Transforms c(r) -> C(k)
        for (int i = 0; i < nodes; i++) {
            double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
            const double *row0 = &K0_fwd[i * nodes];
            const double *row2 = &K2_fwd[i * nodes];
            for (int j = 0; j < nodes; j++) {
                sum0 += row0[j] * c->data[0][j];
                sum1 += row0[j] * c->data[1][j];
                sum2 += row2[j] * c->data[2][j];
            }
            C_k->data[0][i] = sum0;
            C_k->data[1][i] = sum1;
            C_k->data[2][i] = sum2;
        }

        // B. Solve OZ in k-space
        solve_oz_k_space(C_k->data, H_k->data, nodes, rho);

        // C. Inverse Hankel Transforms H(k) -> h(r)
        for (int i = 0; i < nodes; i++) {
            double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
            const double *row0 = &K0_inv[i * nodes];
            const double *row2 = &K2_inv[i * nodes];
            for (int j = 0; j < nodes; j++) {
                sum0 += row0[j] * H_k->data[0][j];
                sum1 += row0[j] * H_k->data[1][j];
                sum2 += row2[j] * H_k->data[2][j];
            }
            h->data[0][i] = sum0;
            h->data[1][i] = sum1;
            h->data[2][i] = sum2;
        }

        // D. Calculate Eta = h - c
        for (int p = 0; p < n_projections; p++) {
            for (int i = 0; i < nodes; i++) {
                eta->data[p][i] = h->data[p][i] - c->data[p][i];
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

        if (iter % 10 == 0 || error <= tolerance) {
            printf("Iter %4d: Error = %.5e\n", iter, error);
        }

        if (error <= tolerance) {
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

        // Apply Anderson Mixing if history >= 2
        int applied_anderson = 0;
        if (hist_count >= 2) {
            int m = hist_count;
            double Amat[6][6] = {0};
            double bvec[6] = {0};
            double xvec[6] = {0};

            // Compute Gram matrix G_ab = <d_a, d_b>
            for (int a = 0; a < m; a++) {
                for (int b = a; b < m; b++) {
                    double dot = 0.0;
                    for (int k_idx = 0; k_idx < total_dim; k_idx++) {
                        dot += d_hist[a][k_idx] * d_hist[b][k_idx];
                    }
                    Amat[a][b] = dot;
                    Amat[b][a] = dot;
                }
                Amat[a][m] = 1.0;
                Amat[m][a] = 1.0;
            }
            Amat[m][m] = 0.0;
            bvec[m] = 1.0;

            if (solve_small_system(m + 1, Amat, bvec, xvec)) {
                // Check stability of weights
                double weight_sum = 0.0;
                int stable = 1;
                for (int j = 0; j < m; j++) {
                    weight_sum += fabs(xvec[j]);
                    if (fabs(xvec[j]) > 50.0) stable = 0;
                }
                if (stable && weight_sum < 100.0) {
                    idx = 0;
                    for (int p = 0; p < n_projections; p++) {
                        for (int i = 0; i < nodes; i++) {
                            double val = 0.0;
                            for (int j = 0; j < m; j++) {
                                val += xvec[j] * (c_hist[j][idx] + alpha * d_hist[j][idx]);
                            }
                            c->data[p][i] = val;
                            idx++;
                        }
                    }
                    applied_anderson = 1;
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

        free(d_curr);
        free(c_curr);
        iter++;
    }

    printf("Iter %4d: Error = %.5e  [DONE]\n", iter, error);
    free_projection_matrix(c_new_mat);

    // Free Anderson history
    for (int m = 0; m < M_DEPTH; m++) {
        free(c_hist[m]);
        free(d_hist[m]);
    }

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
