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
                c[0][i] = h000 - log(g000) + quad_term;
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
 * Evaluates the integro-differential formulation using the exact (000, 110, 112) algebraic reduction.
 */
void closure_RHNC_dipolar(double **c, double **h, double **eta, double *r, int n_points, double beta_mu2, double sigma, double *c_HS, double *h_HS) {
    double dr = r[1] - r[0]; // Assume uniform grid

    // Allocate arrays for the integrand I_k(r)
    double *I0 = malloc(n_points * sizeof(double));
    double *I1 = malloc(n_points * sizeof(double));
    double *I2 = malloc(n_points * sizeof(double));

    for (int i = 0; i < n_points; i++) {
        double ri = r[i];
        if (ri <= sigma) {
            I0[i] = 0.0; I1[i] = 0.0; I2[i] = 0.0;
            continue; // Inside core handled centrally later
        }

        // Compute derivatives using central differences
        int i_prev = (i > 0) ? i - 1 : 0;
        int i_next = (i < n_points - 1) ? i + 1 : n_points - 1;
        double r_prev = r[i_prev], r_next = r[i_next];
        double den = r_next - r_prev;
        if (den < 1e-12) den = 2.0*dr;

        // Delta W = -Delta eta + beta Delta u
        // Delta eta_0 = eta0 - (h_HS - c_HS)
        double eta_HS_prev = h_HS[i_prev] - c_HS[i_prev];
        double eta_HS_next = h_HS[i_next] - c_HS[i_next];
        double dW0_prev = -(eta[0][i_prev] - eta_HS_prev);
        double dW0_next = -(eta[0][i_next] - eta_HS_next);
        double dW0 = (dW0_next - dW0_prev) / den;

        double dW1 = -(eta[1][i_next] - eta[1][i_prev]) / den;

        // u_D = - beta_mu2 / r^3 => beta Delta u_2 = beta_mu2 / r^3
        double u2_prev = beta_mu2 / pow(r_prev, 3.0);
        double u2_next = beta_mu2 / pow(r_next, 3.0);
        double dW2_prev = -eta[2][i_prev] + u2_prev;
        double dW2_next = -eta[2][i_next] + u2_next;
        double dW2 = (dW2_next - dW2_prev) / den;

        // dW_HS = d(ln(g_HS))/dr
        double g_HS_prev = h_HS[i_prev] + 1.0;
        double g_HS_next = h_HS[i_next] + 1.0;
        if (g_HS_prev < 1e-12) g_HS_prev = 1e-12;
        if (g_HS_next < 1e-12) g_HS_next = 1e-12;
        double dW_HS = (log(g_HS_next) - log(g_HS_prev)) / den;

        // Values at current i
        double h0 = h[0][i];
        double h1 = h[1][i];
        double h2 = h[2][i];
        double dh0 = h0 - h_HS[i];

        // Exact RHNC Tensor multiplication rules for l<=2
        I0[i] = dW0 * h0 + 3.0 * dW1 * h1 + (2.0/3.0) * dW2 * h2 - dh0 * dW_HS;
        I1[i] = dW0 * h1 + dW1 * h0 - h1 * dW_HS;
        I2[i] = dW0 * h2 + dW2 * h0 - h2 * dW_HS;
    }

    // Integrate backward: Delta c_k(r) = Integral_r^infty I_k(r') dr'
    double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
    for (int i = n_points - 1; i >= 0; i--) {
        if (r[i] > sigma) {
            if (i < n_points - 1) {
                // Trapezoidal rule
                sum0 += 0.5 * (I0[i] + I0[i+1]) * dr;
                sum1 += 0.5 * (I1[i] + I1[i+1]) * dr;
                sum2 += 0.5 * (I2[i] + I2[i+1]) * dr;
            } else {
                sum0 += I0[i] * dr;
                sum1 += I1[i] * dr;
                sum2 += I2[i] * dr;
            }

            // c = Delta c + c_HS - beta Delta u
            c[0][i] = sum0 + c_HS[i];
            c[1][i] = sum1;
            c[2][i] = sum2 + (beta_mu2 / pow(r[i], 3.0));
        } else {
            // Inside Hard Core exact relation: c = -1 - eta
            c[0][i] = -1.0 - eta[0][i];
            c[1][i] = -eta[1][i];
            c[2][i] = -eta[2][i];
        }
    }

    free(I0);
    free(I1);
    free(I2);
}
