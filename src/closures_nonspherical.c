#include "structures_nonspherical.h"
#include <math.h>
#include <stdio.h>

/**
 * @brief Applies the Mean Spherical Approximation (MSA) closure.
 */
// Updated MSA to accept eta and enforce Core condition
void closure_MSA_dipolar(double **c, double **eta, double *r, int n_points, double beta_mu2, double sigma) {
    for (int i = 0; i < n_points; i++) {
        double ri = r[i];
        if (ri > sigma) {
            c[0][i] = 0.0; 
            c[1][i] = 0.0;
            c[2][i] = beta_mu2 / pow(ri, 3.0);
        } else {
            // Inside Hard Core: h(r) = -1 => c(r) = -1 - eta(r)
            // Ideally we need eta here.
            c[0][i] = -1.0 - eta[0][i];
            c[1][i] = -eta[1][i]; // h110 = 0 inside? usually h(r) = -1 refers to g(r)=0. 
                                  // g = h + 1. So g=0 => h=-1.
                                  // For projections: g000 = h000 + 1. g110 = h110. g112 = h112.
                                  // Condition g(r, omega) = 0 for all omega inside core.
                                  // => h000 = -1, h110 = 0, h112 = 0.
            c[2][i] = -eta[2][i];
        }
    }
}

/**
 * @brief Applies the Linearized HNC (LHNC) closure.
 */
void closure_LHNC_dipolar(double **c, double **h, double **eta, double *r, int n_points, double beta_mu2, double sigma) {
    for (int i = 0; i < n_points; i++) {
        double ri = r[i];
        if (ri > sigma) {
            double h000 = h[0][i];
            double g000 = h000 + 1.0;
            double dipole = beta_mu2 / pow(ri, 3.0);

            // c000
            if (g000 > 1e-12) {
                c[0][i] = h000 - log(g000); 
            } else {
                c[0][i] = -1.0; // Fallback
            }

            // c110
            c[1][i] = h000 * eta[1][i];

            // c112
            c[2][i] = dipole + h000 * (eta[2][i] + dipole);
        } else {
             // Inside Core: Exact relation for HNC is g=0 => c = -1 - eta (same as MSA/PercusYevick approximation inside)
             // or solves for c such that g=0.
             // "In the mean spherical approximation (MSA) the closure relation is ... h(r) = -1 for r < sigma"
             // For HNC, we also enforce h = -1.
             // If we just use the HNC expression c = h - ln(g) + ... it diverges if g->0.
             // So for numerical solvers one usually switches to:
             // if r < sigma: c = -1 - eta (enforces h=-1)
             // else: c = HNC_expression
             
            c[0][i] = -1.0 - eta[0][i];
            c[1][i] = -eta[1][i]; 
            c[2][i] = -eta[2][i];
        }
    }
}

/**
 * @brief Applies the Quadratic HNC (QHNC) closure.
 */
void closure_QHNC_dipolar(double **c, double **h, double **eta, double *r, int n_points, double beta_mu2, double sigma) {
     for (int i = 0; i < n_points; i++) {
        double ri = r[i];
        if (ri > sigma) {
            double h000 = h[0][i];
            double g000 = h000 + 1.0;
            double dipole = beta_mu2 / pow(ri, 3.0);
            
            double term112 = eta[2][i] + dipole;
            double term110 = eta[1][i];
            double quad_term = (term112*term112 + term110*term110) / 6.0;

            if (g000 > 1e-12) {
                c[0][i] = h000 - log(g000) + log(1.0 + quad_term);
            } else {
                c[0][i] = -1.0;
            }

            c[1][i] = h000 * eta[1][i];
            c[2][i] = dipole + h000 * (term112);
        } else {
            // Inside Core
            c[0][i] = -1.0 - eta[0][i];
            c[1][i] = -eta[1][i]; 
            c[2][i] = -eta[2][i];
        }
    }
}

/**
 * @brief Applies the exact Reference Hypernetted-Chain (RHNC) closure for Dipolar Hard Spheres.
 * Evaluates the integro-differential formulation using the exact rotational invariant projections:
 * (000, 110, 112) of Fries & Patey (1985).
 */
void closure_RHNC_dipolar(double **c, double **h, double **eta, double *r, int n_points, 
                          double beta_mu2, double sigma, 
                          const double *c_HS, const double *h_HS, 
                          const double *d_ln_g_HS, const double *d_eta_HS) {
    double dr = r[1] - r[0];

    // Allocate arrays for total radial derivatives D0, D1, D2
    double *D0 = malloc(n_points * sizeof(double));
    double *D1 = malloc(n_points * sizeof(double));
    double *D2 = malloc(n_points * sizeof(double));

    for (int i = 0; i < n_points; i++) {
        if (r[i] <= sigma) {
            D0[i] = 0.0;
            D1[i] = 0.0;
            D2[i] = 0.0;
            continue;
        }

        // Central difference derivatives for eta[p]
        double d_eta0 = (i < n_points - 1 && i > 0) ? (eta[0][i + 1] - eta[0][i - 1]) / (2.0 * dr) : 0.0;
        double d_eta1 = (i < n_points - 1 && i > 0) ? (eta[1][i + 1] - eta[1][i - 1]) / (2.0 * dr) : 0.0;
        double d_eta2 = (i < n_points - 1 && i > 0) ? (eta[2][i + 1] - eta[2][i - 1]) / (2.0 * dr) : 0.0;

        // X' = d(Delta eta - beta Delta u)/dr
        double r4 = pow(r[i], 4.0);
        double X_prime0 = d_eta0 - d_eta_HS[i];
        double X_prime1 = d_eta1;
        double X_prime2 = d_eta2 - (3.0 * beta_mu2 / r4);

        double h0 = h[0][i];
        double h1 = h[1][i];
        double h2 = h[2][i];

        // Total derivatives for r > sigma:
        // D0 = h0 * X'0 + (1/3)*h1*X'1 + (2/3)*h2*X'2 + (h0 - h_HS) * d(ln g_HS)/dr
        // D1 = h0 * X'1 + h1 * X'0 + h1 * d(ln g_HS)/dr
        // D2 = h0 * X'2 + h2 * X'0 + h2 * d(ln g_HS)/dr - 3*beta*mu^2 / r^4
        D0[i] = h0 * X_prime0 + (1.0 / 3.0) * h1 * X_prime1 + (2.0 / 3.0) * h2 * X_prime2 + (h0 - h_HS[i]) * d_ln_g_HS[i];
        D1[i] = h0 * X_prime1 + h1 * X_prime0 + h1 * d_ln_g_HS[i];
        D2[i] = h0 * X_prime2 + h2 * X_prime0 + h2 * d_ln_g_HS[i] - (3.0 * beta_mu2 / r4);
    }

    // Backwards integration from r_max down to sigma:
    // Delta c(r) = - \int_r^{r_max} D(s) ds
    double Delta_c0 = 0.0;
    double Delta_c1 = 0.0;
    double Delta_c2 = 0.0;

    for (int i = n_points - 1; i >= 0; i--) {
        if (r[i] > sigma) {
            if (i < n_points - 1) {
                Delta_c0 -= 0.5 * dr * (D0[i] + D0[i + 1]);
                Delta_c1 -= 0.5 * dr * (D1[i] + D1[i + 1]);
                Delta_c2 -= 0.5 * dr * (D2[i] + D2[i + 1]);
            }
            c[0][i] = Delta_c0 + c_HS[i];
            c[1][i] = Delta_c1;
            c[2][i] = Delta_c2;
        } else {
            // Inside hard core: h(r) = -1 => c = -1 - eta for 000, c = -eta for 110, 112
            c[0][i] = -1.0 - eta[0][i];
            c[1][i] = -eta[1][i];
            c[2][i] = -eta[2][i];
        }
    }

    free(D0);
    free(D1);
    free(D2);
}

