#include "structures_nonspherical.h"
#include "facdes2Y.h"
#include "math_aux.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Forward declarations of closure functions
void closure_MSA_dipolar(double **c, double **eta, double *r, int n_points, double beta_mu2, double sigma);
void closure_LHNC_dipolar(double **c, double **h, double **eta, double *r, int n_points, double beta_mu2, double sigma);
void closure_QHNC_dipolar(double **c, double **h, double **eta, double *r, int n_points, double beta_mu2, double sigma);
void closure_RHNC_dipolar(double **c, double **h, double **eta, double *r, int n_points, double beta_mu2, double sigma, double *c_HS, double *h_HS);

/**
 * @brief Solves the OZ equation in k-space for Dipolar Hard Spheres.
 *
 * Uses the chi-basis decoupling (Blum/Wertheim) which reduces the coupled
 * 2x2 system {110, 112} to two independent scalar equations.
 *
 * From the Patey OZ equations:
 *   eta^000 = rho * C000 * H000    (000 decouples)
 *   eta^0   = (rho/3) * C^0 * H^0  (chi=0 mode, C^0 = C110 + 2*C112)
 *   eta^1   = (rho/3) * C^1 * H^1  (chi=1 mode, C^1 = C110 - C112)
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
        double denom1 = 1.0 - (rho / 3.0) * C1;

        double H0 = (fabs(denom0) > 1e-12) ? C0 / denom0 : 0.0;
        double H1 = (fabs(denom1) > 1e-12) ? C1 / denom1 : 0.0;

        // --- Recover Patey projections ---
        H_k[1][i] = (H0 + 2.0 * H1) / 3.0;  // H110
        H_k[2][i] = (H0 - H1) / 3.0;         // H112
    }
}

/**
 * @brief Computes the exact Percus-Yevick Hard Sphere reference functions.
 * Evaluates the exact PY polynomial for c(r), transforms to C(k), solves OZ for H(k),
 * and transforms back to h(r).
 */
void compute_HS_reference(double *c_HS, double *h_HS, double *r, double *k, 
                          int nodes, double dr, double rho, double sigma) {
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

    // 2. Transform to C(k) (Sine transform)
    double *C_k = malloc(nodes * sizeof(double));
    for (int i = 0; i < nodes; i++) {
        double k_val = k[i];
        double sum = 0.0;
        for (int j = 0; j < nodes; j++) {
            double r_val = r[j];
            double j0 = (fabs(k_val*r_val) < 1e-6) ? 1.0 : sin(k_val*r_val)/(k_val*r_val);
            sum += r_val * r_val * c_HS[j] * j0 * dr;
        }
        C_k[i] = 4.0 * M_PI * sum;
    }

    // 3. Solve OZ for H(k)
    double *H_k = malloc(nodes * sizeof(double));
    for (int i = 0; i < nodes; i++) {
        double denom = 1.0 - rho * C_k[i];
        H_k[i] = (fabs(denom) > 1e-12) ? C_k[i] / denom : 0.0;
    }

    // 4. Transform back to h(r) (Inverse Sine transform)
    double dk = k[1] - k[0];
    for (int i = 0; i < nodes; i++) {
        double r_val = r[i];
        double sum = 0.0;
        for (int j = 0; j < nodes; j++) {
            double k_val = k[j];
            double j0 = (fabs(k_val*r_val) < 1e-6) ? 1.0 : sin(k_val*r_val)/(k_val*r_val);
            sum += k_val * k_val * H_k[j] * j0 * dk;
        }
        h_HS[i] = (1.0 / (2.0 * M_PI * M_PI)) * sum;
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
    printf("Closure: %d (0=MSA, 1=LHNC, 2=QHNC)\n", closureID);
    printf("Dipole Moment (mu): %.4f\n", dipole_moment);
    printf("Density (rho): %.4f\n", rho);

    // 1. Initialize Data Structures
    int n_projections = 3; // 000, 110, 112
    ProjectionMatrix *h = create_projection_matrix(n_projections, nodes);
    ProjectionMatrix *c = create_projection_matrix(n_projections, nodes);
    ProjectionMatrix *eta = create_projection_matrix(n_projections, nodes); // eta = h - c
    
    // Arrays for k-space
    ProjectionMatrix *C_k = create_projection_matrix(n_projections, nodes);
    ProjectionMatrix *H_k = create_projection_matrix(n_projections, nodes);

    // Grid generation
    double dr = rmax / nodes;
    double *r = malloc(nodes * sizeof(double));
    double *k = malloc(nodes * sizeof(double));
    double dk = M_PI / (nodes * dr); // Approximate Nyquist (actually PI/dr?) 
    // Usually dk = pi / (N*dr) or similar. 
    // Standard FFT grid: dk = 2*pi / (N*dr). 
    // Let's stick to what main solver likely uses or standard discrete definition.
    // If r[i] = (i+1)*dr, then k[i] = (i+1)*dk.
    dk = M_PI / (nodes * dr);

    for(int i=0; i<nodes; i++) {
        r[i] = (i+1) * dr;
        k[i] = (i+1) * dk;
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
        compute_HS_reference(c_HS, h_HS, r, k, nodes, dr, rho, sigma);
    }

    // 2. Initialization 
    closure_MSA_dipolar(c->data, eta->data, r, nodes, beta_mu2, sigma);

    // 3. Iteration Loop
    int max_iter = 2000;
    double tolerance = 1e-6;
    double error = 1.0;
    int iter = 0;

    printf("Starting Iteration...\n");
    while (iter < max_iter && error > tolerance) {
        
        // A. Transforms c(r) -> C(k)
        // 000: Sine Transform (j0)
        // We use FFT wrapper from math_aux which does sine transform
        // Note: FFT in math_aux acts IN PLACE? Uses 'sinft'. 
        // Copies to C_k first.
        for(int i=0; i<nodes; i++) C_k->data[0][i] = c->data[0][i] * r[i]; // r*c(r) for sine transform?
        // Check FFT implementation. math_aux FFT assumes input is f(r). 
        // Actually math_aux FFT: "inputData[i] = (i * inputData[i]);" inside? No, it does sinft directly.
        // Let's use generic sinft calls if possible, or copy FFT logic. 
        // math_aux's FFT function: "inputData[i] = (i * inputData[i])" ?? That looks suspicious for internal rescaling.
        // Let's use our new HT2 logic for 112, and maybe just call sine transform for 000/110.
        // For 000/110 (l=0) -> 4*pi/k * Integral(r * c(r) * sin(kr))
        // So we transform f(r) = r*c(r) using sine transform.
        
        // Manual Sine Transform for 000 and 110 using simple direct sum for consistency with HT2 approx for now
        // or re-use HT2 with j0? j0(x) = sin(x)/x.
        // Integral r^2 c(r) j0(kr) dr = Integral r^2 c(r) sin(kr)/(kr) dr = 1/k Integral r c(r) sin(kr) dr.
        // This is Sine Transform of r*c(r) divided by k.
        
        // Using HT2_Direct logic but with j0 for 000 and 110
        // NOTE: Efficient FFT not used here yet, keeping O(N^2) for correctness first.
        // We need a HT0_Direct equivalent.
        // Let's just use HT2_Direct but swap j2 for j0... 
        // Hardcoding loop here for now for 000/110
         for (int i = 0; i < nodes; i++) {
            double k_val = k[i];
            double sum0 = 0.0;
            double sum1 = 0.0;
            for (int j = 0; j < nodes; j++) {
                double r_val = r[j];
                double j0 = (fabs(k_val*r_val) < 1e-6) ? 1.0 : sin(k_val*r_val)/(k_val*r_val);
                double term0 = r_val * r_val * c->data[0][j] * j0;
                double term1 = r_val * r_val * c->data[1][j] * j0;
                
                double w = dr; // simple weight
                sum0 += term0 * w;
                sum1 += term1 * w;
            }
            C_k->data[0][i] = 4.0 * M_PI * sum0;
            C_k->data[1][i] = 4.0 * M_PI * sum1;
        }

        // 112: Hankel Transform order 2
        HT2_Direct(c->data[2], C_k->data[2], r, k, nodes);

        // B. Solve OZ in k-space
        solve_oz_k_space(C_k->data, H_k->data, nodes, rho);

        // C. Transforms H(k) -> h(r)
        // 000/110: Inverse Sine / j0
        // h(r) = 1/(2pi^2) Integral k^2 H(k) j0(kr) dk
        for (int i = 0; i < nodes; i++) {
            double r_val = r[i];
            double sum0 = 0.0;
            double sum1 = 0.0;
            for (int j = 0; j < nodes; j++) {
                double k_val = k[j];
                double j0 = (fabs(k_val*r_val) < 1e-6) ? 1.0 : sin(k_val*r_val)/(k_val*r_val);
                double term0 = k_val * k_val * H_k->data[0][j] * j0;
                double term1 = k_val * k_val * H_k->data[1][j] * j0;
                
                double w = dk;
                sum0 += term0 * w;
                sum1 += term1 * w;
            }
            h->data[0][i] = sum0 / (2.0 * M_PI * M_PI);
            h->data[1][i] = sum1 / (2.0 * M_PI * M_PI);
        }

        // 112: Inverse HT2
        IHT2_Direct(H_k->data[2], h->data[2], r, k, nodes);

        // D. Calculate Eta = h - c
        for(int p=0; p<n_projections; p++)
            for(int i=0; i<nodes; i++)
                eta->data[p][i] = h->data[p][i] - c->data[p][i];

        // E. Closure: compute c_new from h and eta
        // Allocate a temporary ProjectionMatrix to hold the closure output
        ProjectionMatrix *c_new_mat = create_projection_matrix(n_projections, nodes);
        for(int p=0; p<n_projections; p++)
            for(int i=0; i<nodes; i++)
                c_new_mat->data[p][i] = c->data[p][i];

        // E. Apply Closure Extension c(r)
        if (closureID == 0) {
            closure_MSA_dipolar(c_new_mat->data, eta->data, r, nodes, beta_mu2, sigma);
        } else if (closureID == 1) {
            closure_LHNC_dipolar(c_new_mat->data, h->data, eta->data, r, nodes, beta_mu2, sigma);
        } else if (closureID == 2) {
            closure_QHNC_dipolar(c_new_mat->data, h->data, eta->data, r, nodes, beta_mu2, sigma);
        } else if (closureID == 3) {
            closure_RHNC_dipolar(c_new_mat->data, h->data, eta->data, r, nodes, beta_mu2, sigma, c_HS, h_HS);
        } 

        // F. Compute L2 error and apply Picard mixing
        // c_mixed = (1 - alpha) * c_old + alpha * c_new
        double alpha = 0.3;   // Damping factor (0 < alpha <= 1)
        error = 0.0;
        for(int p=0; p<n_projections; p++) {
            for(int i=0; i<nodes; i++) {
                double diff = c_new_mat->data[p][i] - c->data[p][i];
                error += diff * diff;
                c->data[p][i] += alpha * diff;  // Picard update
            }
        }
        error = sqrt(error / (n_projections * nodes));
        free_projection_matrix(c_new_mat);

        if (iter % 50 == 0)
            printf("Iter %4d: Error = %.5e\n", iter, error);
        iter++;
    }
    printf("Iter %4d: Error = %.5e  [DONE]\n", iter-1, error);

    // 4. Output Results
    FILE *fp = fopen("output/output_dipolar.dat", "w");
    if(fp) {
        fprintf(fp, "# r h000 h110 h112 c000 c110 c112\n");
        for(int i=0; i<nodes; i++) {
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
    if(fp_k) {
        fprintf(fp_k, "# k H000 H110 H112 C000 C110 C112\n");
        for(int i=0; i<nodes; i++) {
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
    // Patey:
    //   S000 = 1 + rho * H000   (equivalently: 1 / (1 - rho*C000))
    //   S110 = rho * H110        (or via chi inversion)
    //   S112 = rho * H112        (or via chi inversion)
    //
    // Chi representation:
    //   C0 = C110 + 2*C112,   S0 = 1 / (1 - rho/3 * C0)
    //   C1 = C110 -   C112,   S1 = 1 / (1 + rho/3 * C1)
    //
    // Recovery (chi -> Patey):
    //   S110 = (S0 + 2*S1)/3 - 1
    //   S112 = (S0 - S1)/3
    // ---------------------------------------------------------------
    FILE *fp_sk = fopen("output/output_dipolar_sk.dat", "w");
    if(fp_sk) {
        fprintf(fp_sk, "# k  S000  S110_Patey  S112_Patey  S0_chi  S1_chi\n");
        for(int i = 0; i < nodes; i++) {
            double C000 = C_k->data[0][i];
            double C110 = C_k->data[1][i];
            double C112 = C_k->data[2][i];

            // --- Patey (direct from H_k) ---
            double S000 = 1.0 + rho * H_k->data[0][i];
            double S110 = rho * H_k->data[1][i];
            double S112 = rho * H_k->data[2][i];

            // --- Chi representation (from C_k) ---
            double C0 = C110 + 2.0 * C112;  // chi=0 mode
            double C1 = C110 - C112;         // chi=1 mode

            double denom0 = 1.0 - (rho / 3.0) * C0;
            double denom1 = 1.0 + (rho / 3.0) * C1;

            double S0 = (fabs(denom0) > 1e-12) ? 1.0 / denom0 : 1e12;
            double S1 = (fabs(denom1) > 1e-12) ? 1.0 / denom1 : 1e12;

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
    if (c_HS) free(c_HS);
    if (h_HS) free(h_HS);
}

