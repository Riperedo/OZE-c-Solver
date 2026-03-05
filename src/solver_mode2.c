#include "structures_nonspherical.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// ----------------------------------------------------
// Explicit Discrete Sine/Cosine Hankel Transforms
// ----------------------------------------------------

// Fast approximation using explicit Bessel expansions to avoid orthogonality drift in Picard.
static double get_jl_kr(int l, double kr) {
    if (kr < 1e-12) {
        return (l == 0) ? 1.0 : 0.0;
    }
    double kr2 = kr * kr;
    double kr3 = kr2 * kr;
    double kr4 = kr3 * kr;
    double kr5 = kr4 * kr;
    
    double sk = sin(kr);
    double ck = cos(kr);
    
    switch(l) {
        case 0: return sk / kr;
        case 1: return (sk / kr2) - (ck / kr);
        case 2: return ((3.0 / kr3) - (1.0 / kr)) * sk - (3.0 / kr2) * ck;
        case 3: return ((15.0 / kr4) - (6.0 / kr2)) * sk + ((1.0 / kr) - (15.0 / kr3)) * ck;
        case 4: return ((105.0 / kr5) - (45.0 / kr3) + (1.0 / kr)) * sk + ((10.0 / kr2) - (105.0 / kr4)) * ck;
        default: return 0.0;
    }
}

// Modes derived from m,n <= 2 and m+n+l even
// Total 14 modes for symmetric non-spherical potentials like dipoles
static int mode_l[14] = {
    0, 1, 2, 1, 0, 2, 1, 3, 2, 1, 3, 0, 2, 4
};
// Correct list ordered identically to Python generator:
// 0: 000(l=0), 1: 011(1), 2: 022(2), 3: 101(1), 4: 110(0), 5: 112(2),
// 6: 121(1), 7: 123(3), 8: 202(2), 9: 211(1), 10: 213(3),
// 11: 220(0), 12: 222(2), 13: 224(4)

static int get_mode_l(int index) {
    return mode_l[index];
}

void solve_oz_k_space_mode2(ProjectionMatrix *C_mat, ProjectionMatrix *H_mat, int nodes, double rho);

void closure_MSA_mode2(double **c, double **eta, double *r, int n_points, double beta_mu2, double sigma, int n_projections) {
    for (int i = 0; i < n_points; i++) {
        if (r[i] > sigma) {
            for(int p=0; p<n_projections; p++) {
                c[p][i] = 0.0;
            }
            c[5][i] = beta_mu2 / pow(r[i], 3.0); // Index 5 is 112
        } else {
            c[0][i] = -1.0 - eta[0][i];
            for(int p=1; p<n_projections; p++) {
                c[p][i] = -eta[p][i];
            }
        }
    }
}

void closure_LHNC_mode2(double **c, double **h, double **eta, double *r, int n_points, double beta_mu2, double sigma, int n_projections) {
    for (int i = 0; i < n_points; i++) {
        if (r[i] > sigma) {
            double h000 = h[0][i];
            double g000 = h000 + 1.0;
            if (g000 < 1e-12) g000 = 1e-12;
            c[0][i] = h000 - log(g000);
            
            for(int p=1; p<n_projections; p++) {
                c[p][i] = h000 * eta[p][i];
            }
            c[5][i] += beta_mu2 / pow(r[i], 3.0) + h000 * (beta_mu2 / pow(r[i], 3.0));
        } else {
            c[0][i] = -1.0 - eta[0][i];
            for(int p=1; p<n_projections; p++) {
                c[p][i] = -eta[p][i];
            }
        }
    }
}

// Quick integration using generic spherical bessel
void solver_mode2_core(int closureID, double temp, double rho, double dipole_moment, 
                   int nodes, double rmax, const char *output_dir) {
    
    printf("Initializing Extended Mode 2 Solver (Potential 15, m,n<=2, parity even)...\n");
    int n_projections = 14; 
    ProjectionMatrix *h = create_projection_matrix(n_projections, nodes);
    ProjectionMatrix *c = create_projection_matrix(n_projections, nodes);
    ProjectionMatrix *eta = create_projection_matrix(n_projections, nodes); 
    ProjectionMatrix *C_k = create_projection_matrix(n_projections, nodes);
    ProjectionMatrix *H_k = create_projection_matrix(n_projections, nodes);

    double dr = rmax / nodes;
    double *r = malloc(nodes * sizeof(double));
    double *k = malloc(nodes * sizeof(double));
    double dk = M_PI / (nodes * dr); 

    for(int i=0; i<nodes; i++) {
        r[i] = (i+1) * dr;
        k[i] = (i+1) * dk;
    }

    double beta = 1.0 / temp;
    double beta_mu2 = beta * dipole_moment * dipole_moment;
    double sigma = 1.0;

    closure_MSA_mode2(c->data, eta->data, r, nodes, beta_mu2, sigma, n_projections);

    int max_iter = 2000;
    double tolerance = 1e-6;
    double error = 1.0;
    int iter = 0;

    while (iter < max_iter && error > tolerance) {
        // Forward Hankel Transform O(N^2)
        for (int p = 0; p < n_projections; p++) {
            int l = get_mode_l(p);
            for (int i = 0; i < nodes; i++) {
                double k_val = k[i];
                double sum = 0.0;
                for (int j = 0; j < nodes; j++) {
                    double r_val = r[j];
                    double arg = k_val * r_val;
                    double bessel = (arg < 1e-6 && l > 0) ? 0.0 : ((arg < 1e-6 && l == 0) ? 1.0 : get_jl_kr(l, arg));
                    sum += r_val * r_val * c->data[p][j] * bessel * dr;
                }
                C_k->data[p][i] = 4.0 * M_PI * sum;
            }
        }

        solve_oz_k_space_mode2(C_k, H_k, nodes, rho);

        // Inverse Hankel Transform O(N^2)
        for (int p = 0; p < n_projections; p++) {
            int l = get_mode_l(p);
            for (int i = 0; i < nodes; i++) {
                double r_val = r[i];
                double sum = 0.0;
                for (int j = 0; j < nodes; j++) {
                    double k_val = k[j];
                    double arg = k_val * r_val;
                    double bessel = (arg < 1e-6 && l > 0) ? 0.0 : ((arg < 1e-6 && l == 0) ? 1.0 : get_jl_kr(l, arg));
                    sum += k_val * k_val * H_k->data[p][j] * bessel * dk;
                }
                h->data[p][i] = (1.0 / (2.0 * M_PI * M_PI)) * sum;
                eta->data[p][i] = h->data[p][i] - c->data[p][i];
            }
        }

        ProjectionMatrix *c_new = create_projection_matrix(n_projections, nodes);
        for(int p=0; p<n_projections; p++) {
            for(int i=0; i<nodes; i++) {
                c_new->data[p][i] = c->data[p][i];
            }
        }

        if (closureID == 0) closure_MSA_mode2(c_new->data, eta->data, r, nodes, beta_mu2, sigma, n_projections);
        else if (closureID == 1) closure_LHNC_mode2(c_new->data, h->data, eta->data, r, nodes, beta_mu2, sigma, n_projections);
        else closure_LHNC_mode2(c_new->data, h->data, eta->data, r, nodes, beta_mu2, sigma, n_projections); // Fallback

        double alpha = 0.3;
        error = 0.0;
        for(int p=0; p<n_projections; p++) {
            for(int i=0; i<nodes; i++) {
                double diff = c_new->data[p][i] - c->data[p][i];
                error += diff * diff;
                c->data[p][i] += alpha * diff;
            }
        }
        error = sqrt(error / (n_projections * nodes));
        free_projection_matrix(c_new);

        if (iter % 50 == 0) printf("Iter %4d: Error = %.5e\n", iter, error);
        iter++;
    }
    printf("Finished Mode 2 Solver in %d iter. Error = %.5e\n", iter-1, error);
    
    // Save output...
    char filename[256];
    snprintf(filename, sizeof(filename), "output/output_mode15.dat");
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "# r h000 h110 h112 h220 h222 c000 ... \n");
    for (int i=0; i<nodes; i++) {
        fprintf(fp, "%.5e", r[i]);
        for(int p=0; p<n_projections; p++) fprintf(fp, " %.5e", h->data[p][i]);
        fprintf(fp, "\n");
    }
    fclose(fp);

    free_projection_matrix(h); free_projection_matrix(c); free_projection_matrix(eta);
    free_projection_matrix(C_k); free_projection_matrix(H_k);
    free(r); free(k);
}

void solve_oz_k_space_mode2(ProjectionMatrix *C_mat, ProjectionMatrix *H_mat, int nodes, double rho) {
    double **C_k = C_mat->data;
    double **H_k = H_mat->data;
    for (int i = 0; i < nodes; i++) {
        double C_0_Blum = C_k[0][i] / (1.0);
        double C_1_Blum = C_k[1][i] / (-0.9999999999999999);
        double C_2_Blum = C_k[2][i] / (1.0);
        double C_3_Blum = C_k[3][i] / (-0.9999999999999999);
        double C_4_Blum = C_k[4][i] / (-1.732050807568877);
        double C_5_Blum = C_k[5][i] / (1.095445115010332);
        double C_6_Blum = C_k[6][i] / (1.414213562373095);
        double C_7_Blum = C_k[7][i] / (-1.1338934190276817);
        double C_8_Blum = C_k[8][i] / (1.0);
        double C_9_Blum = C_k[9][i] / (1.414213562373095);
        double C_10_Blum = C_k[10][i] / (-1.1338934190276817);
        double C_11_Blum = C_k[11][i] / (2.2360679774997902);
        double C_12_Blum = C_k[12][i] / (-1.1952286093343938);
        double C_13_Blum = C_k[13][i] / (1.1952286093343938);
        double C_p0_00 = (1.0) * C_0_Blum;
        double C_p0_01 = (-0.5773502691896257) * C_1_Blum;
        double C_p0_02 = (0.4472135954999579) * C_2_Blum;
        double C_p0_10 = (-0.5773502691896257) * C_3_Blum;
        double C_m1_11 = (0.5773502691896257) * C_4_Blum + (0.18257418583505536) * C_5_Blum;
        double C_p0_11 = (-0.5773502691896257) * C_4_Blum + (0.3651483716701107) * C_5_Blum;
        double C_p1_11 = (0.5773502691896257) * C_4_Blum + (0.18257418583505536) * C_5_Blum;
        double C_m1_12 = (-0.31622776601683794) * C_6_Blum + (-0.1690308509457033) * C_7_Blum;
        double C_p0_12 = (0.3651483716701107) * C_6_Blum + (-0.29277002188455997) * C_7_Blum;
        double C_p1_12 = (-0.31622776601683794) * C_6_Blum + (-0.1690308509457033) * C_7_Blum;
        double C_p0_20 = (0.4472135954999579) * C_8_Blum;
        double C_m1_21 = (-0.31622776601683794) * C_9_Blum + (-0.1690308509457033) * C_10_Blum;
        double C_p0_21 = (0.3651483716701107) * C_9_Blum + (-0.29277002188455997) * C_10_Blum;
        double C_p1_21 = (-0.31622776601683794) * C_9_Blum + (-0.1690308509457033) * C_10_Blum;
        double C_m2_22 = (0.4472135954999579) * C_11_Blum + (0.23904572186687872) * C_12_Blum + (0.03984095364447979) * C_13_Blum;
        double C_m1_22 = (-0.4472135954999579) * C_11_Blum + (0.11952286093343936) * C_12_Blum + (0.15936381457791915) * C_13_Blum;
        double C_p0_22 = (0.4472135954999579) * C_11_Blum + (-0.23904572186687872) * C_12_Blum + (0.23904572186687872) * C_13_Blum;
        double C_p1_22 = (-0.4472135954999579) * C_11_Blum + (0.11952286093343936) * C_12_Blum + (0.15936381457791915) * C_13_Blum;
        double C_p2_22 = (0.4472135954999579) * C_11_Blum + (0.23904572186687872) * C_12_Blum + (0.03984095364447979) * C_13_Blum;

        // --- Solvers for each Chi mode ---
        double I_m_C_p0_00 = 1.0 - (1.0 * rho) * C_p0_00;
        double I_m_C_p0_01 = 0.0 - (1.0 * rho) * C_p0_01;
        double I_m_C_p0_02 = 0.0 - (1.0 * rho) * C_p0_02;
        double I_m_C_p0_10 = 0.0 - (1.0 * rho) * C_p0_10;
        double I_m_C_p0_11 = 1.0 - (1.0 * rho) * C_p0_11;
        double I_m_C_p0_12 = 0.0 - (1.0 * rho) * C_p0_12;
        double I_m_C_p0_20 = 0.0 - (1.0 * rho) * C_p0_20;
        double I_m_C_p0_21 = 0.0 - (1.0 * rho) * C_p0_21;
        double I_m_C_p0_22 = 1.0 - (1.0 * rho) * C_p0_22;

        double det_p0 = I_m_C_p0_00 * (I_m_C_p0_11 * I_m_C_p0_22 - I_m_C_p0_12 * I_m_C_p0_21)
                              - I_m_C_p0_01 * (I_m_C_p0_10 * I_m_C_p0_22 - I_m_C_p0_12 * I_m_C_p0_20)
                              + I_m_C_p0_02 * (I_m_C_p0_10 * I_m_C_p0_21 - I_m_C_p0_11 * I_m_C_p0_20);
        double inv_det_p0 = (fabs(det_p0) > 1e-12) ? 1.0 / det_p0 : 0.0;

        double M_p0_00 = (I_m_C_p0_11 * I_m_C_p0_22 - I_m_C_p0_12 * I_m_C_p0_21) * inv_det_p0;
        double M_p0_01 = (I_m_C_p0_02 * I_m_C_p0_21 - I_m_C_p0_01 * I_m_C_p0_22) * inv_det_p0;
        double M_p0_02 = (I_m_C_p0_01 * I_m_C_p0_12 - I_m_C_p0_02 * I_m_C_p0_11) * inv_det_p0;
        double M_p0_10 = (I_m_C_p0_12 * I_m_C_p0_20 - I_m_C_p0_10 * I_m_C_p0_22) * inv_det_p0;
        double M_p0_11 = (I_m_C_p0_00 * I_m_C_p0_22 - I_m_C_p0_02 * I_m_C_p0_20) * inv_det_p0;
        double M_p0_12 = (I_m_C_p0_02 * I_m_C_p0_10 - I_m_C_p0_00 * I_m_C_p0_12) * inv_det_p0;
        double M_p0_20 = (I_m_C_p0_10 * I_m_C_p0_21 - I_m_C_p0_11 * I_m_C_p0_20) * inv_det_p0;
        double M_p0_21 = (I_m_C_p0_01 * I_m_C_p0_20 - I_m_C_p0_00 * I_m_C_p0_21) * inv_det_p0;
        double M_p0_22 = (I_m_C_p0_00 * I_m_C_p0_11 - I_m_C_p0_01 * I_m_C_p0_10) * inv_det_p0;
        double H_p0_00 = M_p0_00 * C_p0_00 + M_p0_01 * C_p0_10 + M_p0_02 * C_p0_20;
        double H_p0_01 = M_p0_00 * C_p0_01 + M_p0_01 * C_p0_11 + M_p0_02 * C_p0_21;
        double H_p0_02 = M_p0_00 * C_p0_02 + M_p0_01 * C_p0_12 + M_p0_02 * C_p0_22;
        double H_p0_10 = M_p0_10 * C_p0_00 + M_p0_11 * C_p0_10 + M_p0_12 * C_p0_20;
        double H_p0_11 = M_p0_10 * C_p0_01 + M_p0_11 * C_p0_11 + M_p0_12 * C_p0_21;
        double H_p0_12 = M_p0_10 * C_p0_02 + M_p0_11 * C_p0_12 + M_p0_12 * C_p0_22;
        double H_p0_20 = M_p0_20 * C_p0_00 + M_p0_21 * C_p0_10 + M_p0_22 * C_p0_20;
        double H_p0_21 = M_p0_20 * C_p0_01 + M_p0_21 * C_p0_11 + M_p0_22 * C_p0_21;
        double H_p0_22 = M_p0_20 * C_p0_02 + M_p0_21 * C_p0_12 + M_p0_22 * C_p0_22;
        double I_m_C_p1_11 = 1.0 - (-1.0 * rho) * C_p1_11;
        double I_m_C_p1_12 = 0.0 - (-1.0 * rho) * C_p1_12;
        double I_m_C_p1_21 = 0.0 - (-1.0 * rho) * C_p1_21;
        double I_m_C_p1_22 = 1.0 - (-1.0 * rho) * C_p1_22;

        double det_p1 = I_m_C_p1_11 * I_m_C_p1_22 - I_m_C_p1_12 * I_m_C_p1_21;
        double inv_det_p1 = (fabs(det_p1) > 1e-12) ? 1.0 / det_p1 : 0.0;
        double M_p1_11 = I_m_C_p1_22 * inv_det_p1;
        double M_p1_12 = -I_m_C_p1_12 * inv_det_p1;
        double M_p1_21 = -I_m_C_p1_21 * inv_det_p1;
        double M_p1_22 = I_m_C_p1_11 * inv_det_p1;
        double H_p1_11 = M_p1_11 * C_p1_11 + M_p1_12 * C_p1_21;
        double H_p1_12 = M_p1_11 * C_p1_12 + M_p1_12 * C_p1_22;
        double H_p1_21 = M_p1_21 * C_p1_11 + M_p1_22 * C_p1_21;
        double H_p1_22 = M_p1_21 * C_p1_12 + M_p1_22 * C_p1_22;
        double I_m_C_m1_11 = 1.0 - (-1.0 * rho) * C_m1_11;
        double I_m_C_m1_12 = 0.0 - (-1.0 * rho) * C_m1_12;
        double I_m_C_m1_21 = 0.0 - (-1.0 * rho) * C_m1_21;
        double I_m_C_m1_22 = 1.0 - (-1.0 * rho) * C_m1_22;

        double det_m1 = I_m_C_m1_11 * I_m_C_m1_22 - I_m_C_m1_12 * I_m_C_m1_21;
        double inv_det_m1 = (fabs(det_m1) > 1e-12) ? 1.0 / det_m1 : 0.0;
        double M_m1_11 = I_m_C_m1_22 * inv_det_m1;
        double M_m1_12 = -I_m_C_m1_12 * inv_det_m1;
        double M_m1_21 = -I_m_C_m1_21 * inv_det_m1;
        double M_m1_22 = I_m_C_m1_11 * inv_det_m1;
        double H_m1_11 = M_m1_11 * C_m1_11 + M_m1_12 * C_m1_21;
        double H_m1_12 = M_m1_11 * C_m1_12 + M_m1_12 * C_m1_22;
        double H_m1_21 = M_m1_21 * C_m1_11 + M_m1_22 * C_m1_21;
        double H_m1_22 = M_m1_21 * C_m1_12 + M_m1_22 * C_m1_22;
        double denom_p2 = 1.0 - (1.0 * rho) * C_p2_22;
        double H_p2_22 = (fabs(denom_p2) > 1e-12) ? C_p2_22 / denom_p2 : 0.0;
        double denom_m2 = 1.0 - (1.0 * rho) * C_m2_22;
        double H_m2_22 = (fabs(denom_m2) > 1e-12) ? C_m2_22 / denom_m2 : 0.0;

        // --- Inverse Chi Transforms ---
        double H_0_Blum = (1.0) * H_p0_00;
        double H_1_Blum = (-1.7320508075688772) * H_p0_01;
        double H_2_Blum = (2.23606797749979) * H_p0_02;
        double H_3_Blum = (-1.7320508075688772) * H_p0_10;
        double H_4_Blum = (0.5773502691896257) * H_m1_11 + (-0.5773502691896257) * H_p0_11 + (0.5773502691896257) * H_p1_11;
        double H_5_Blum = (0.9128709291752768) * H_m1_11 + (1.8257418583505536) * H_p0_11 + (0.9128709291752768) * H_p1_11;
        double H_6_Blum = (-0.9486832980505138) * H_m1_12 + (1.0954451150103321) * H_p0_12 + (-0.9486832980505138) * H_p1_12;
        double H_7_Blum = (-1.1832159566199232) * H_m1_12 + (-2.04939015319192) * H_p0_12 + (-1.1832159566199232) * H_p1_12;
        double H_8_Blum = (2.23606797749979) * H_p0_20;
        double H_9_Blum = (-0.9486832980505138) * H_m1_21 + (1.0954451150103321) * H_p0_21 + (-0.9486832980505138) * H_p1_21;
        double H_10_Blum = (-1.1832159566199232) * H_m1_21 + (-2.04939015319192) * H_p0_21 + (-1.1832159566199232) * H_p1_21;
        double H_11_Blum = (0.4472135954999579) * H_m2_22 + (-0.4472135954999579) * H_m1_22 + (0.4472135954999579) * H_p0_22 + (-0.4472135954999579) * H_p1_22 + (0.4472135954999579) * H_p2_22;
        double H_12_Blum = (1.1952286093343936) * H_m2_22 + (0.5976143046671968) * H_m1_22 + (-1.1952286093343936) * H_p0_22 + (0.5976143046671968) * H_p1_22 + (1.1952286093343936) * H_p2_22;
        double H_13_Blum = (0.3585685828003181) * H_m2_22 + (1.4342743312012725) * H_m1_22 + (2.1514114968019085) * H_p0_22 + (1.4342743312012725) * H_p1_22 + (0.3585685828003181) * H_p2_22;
        H_k[0][i] = H_0_Blum * (1.0);
        H_k[1][i] = H_1_Blum * (-0.9999999999999999);
        H_k[2][i] = H_2_Blum * (1.0);
        H_k[3][i] = H_3_Blum * (-0.9999999999999999);
        H_k[4][i] = H_4_Blum * (-1.732050807568877);
        H_k[5][i] = H_5_Blum * (1.095445115010332);
        H_k[6][i] = H_6_Blum * (1.414213562373095);
        H_k[7][i] = H_7_Blum * (-1.1338934190276817);
        H_k[8][i] = H_8_Blum * (1.0);
        H_k[9][i] = H_9_Blum * (1.414213562373095);
        H_k[10][i] = H_10_Blum * (-1.1338934190276817);
        H_k[11][i] = H_11_Blum * (2.2360679774997902);
        H_k[12][i] = H_12_Blum * (-1.1952286093343938);
        H_k[13][i] = H_13_Blum * (1.1952286093343938);
    }
}
