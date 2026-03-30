/**
 * @file structures.c
 * @brief Implementation of the Ornstein-Zernike equation solver and structure factor calculations.
 *
 * This file contains the core logic for initializing potentials, solving the OZ equation
 * using Ng's method, enforcing closure relations (PY, HNC, RY), and calculating
 * thermodynamic properties.
 */

#include "structures.h"
#include "math_aux.h"

/**
 * @brief Initializes the system parameters and interaction potentials.
 *
 * @param fv Volume fraction.
 * @param xnu Potential parameter (e.g., for LJT).
 * @param especie1 Properties of species 1.
 * @param especie2 Properties of species 2.
 * @param rmax Maximum radial distance.
 * @param potentialID ID of the interaction potential to use.
 */
void input(double fv, double xnu, species especie1, species especie2, double rmax, int potentialID) {

    int i;
    double dq;

    rho = (6.0 / M_PI) * fv;

    printf("\n------------------------------\n");
    printf("VOLUME FRACTION = %lf\n\n", fv);

    double sigma1 = especie1.diameter;
    double sigma2 = especie2.diameter;

    sigmaVec[0] = sigma1;
    sigmaVec[1] = (sigma1 + sigma2) / 2.0;
    sigmaVec[2] = sigma2;

    dq = M_PI / rmax;
    dr = rmax / ((double) nrows);

    for (i = 0; i < nrows; i++) {
        r[i] = i * dr;
        q[i] = i * dq;
    }

    POT(especie1, especie2, potentialID, xnu);
}

/**
 * @brief Calculates the interaction potential U(r) and its derivative Up(r).
 *
 * Initializes the potential arrays `U` and `Up` based on the selected `potentialID`.
 *
 * @param especie1 Properties of species 1.
 * @param especie2 Properties of species 2.
 * @param potentialID ID of the potential.
 * @param xnu Potential parameter.
 */
void POT(species especie1, species especie2, int potentialID, double xnu) {

    int i, k;
    double dmed, rlamb;
    double arg1, arg2, arg3, arg4;
    double *Ua, *Ur, *E, *E2, *z, *z2;

    // Allocate memory for potential calculation arrays
    Ua = malloc(nrows*ncols * sizeof(double));
    Ur = malloc(nrows*ncols * sizeof(double));
    E = malloc(ncols * sizeof(double));
    E2 = malloc(ncols * sizeof(double));
    z = malloc(ncols * sizeof(double));
    z2 = malloc(ncols * sizeof(double));

    if (Ua == NULL || Ur == NULL || E == NULL || E2 == NULL || z == NULL || z2 == NULL ) {
        printf("Memory allocation failed in POT.\n");
        return;
    }

    // Mean distance between particles
    dmed = pow(rho, -1.0/3.0);

    switch(potentialID){
        
        case 1:
            // INVERSE POWER LAW: U = T* (sigma/r)^(lambda)
            E[0] = 1.0 / especie1.temperature;
            E[2] = 1.0 / especie2.temperature;
            E[1] = sqrt(E[0] * E[2]);

            z[0] = especie1.lambda;
            z[2] = especie2.lambda;
            z[1] = sqrt(z[0] * z[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < (sigmaVec[k] / 2.0)) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        U[i*ncols + k] = E[k] * pow(sigmaVec[k]/r[i], z[k]);
                        Up[i*ncols + k] = U[i*ncols + k] * (z[k]); // -f(r)*r
                    }
                }
            }

            printf("POTENTIAL:   INVERSE POWER LAW\n\n");
            printf("TEMPERATURE:  %.3lf\n", E[0]);
            printf("POWER:        %.3lf\n", z[0]);
            printf("DIAMETER:     %.3lf\n", sigmaVec[0]);
            printf("------------------------------\n");
            break;
        
        case 2: 
            // TRUNCATED LENNARD-JONES (Repulsive only)
            rlamb = 6.0;

            E[0] = 1.0 / especie1.temperature;
            E[2] = 1.0 / especie2.temperature;
            E[1] = sqrt(E[0] * E[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    arg4 = sigmaVec[k] * pow(2.0, 1.0/rlamb);
                    if ((r[i] < (sigmaVec[k] / 2.0)) || (r[i] > arg4)) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        arg1 = pow(sigmaVec[k] / r[i], rlamb);
                        arg2 = arg1 * arg1;
                        arg3 = 1.0/4.0;
                        U[i*ncols + k] = 4.0 * E[k] * (arg2 - arg1 + arg3);
                        Up[i*ncols + k] = 4.0 * E[k] * rlamb * (2.0*arg2 - arg1);
                    }
                }
            }

            printf("POTENTIAL: TRUNCATED LENNARD-JONES 6-12\n\n");
            printf("TEMPERATURE:  %.3lf\n", E[0]);
            printf("DIAMETER:     %.3lf\n", sigmaVec[0]);
            printf("------------------------------\n");
            break;

        case 3: 
            // TRUNCATED LENNARD-JONES (Minimum at r = sigma)
            E[0] = 1.0 / especie1.temperature;
            E[2] = 1.0 / especie2.temperature;
            E[1] = sqrt(E[0] * E[2]);
            
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if ((r[i] < (sigmaVec[k] / 2.0)) || (r[i] > sigmaVec[k])) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        arg1 = pow(sigmaVec[k] / r[i], xnu);
                        arg2 = arg1 * arg1;
                        arg3 = 1.0;
                        U[i*ncols + k] = E[k] * (arg2 - 2.0*arg1 + arg3);
                        Up[i*ncols + k] = E[k] * 2.0*xnu * (arg2 - arg1);
                    }
                }
            }

            printf("POTENTIAL: TRUNCATED LENNARD-JONES %.1lf-%.1lf\n\n", xnu, 2.0*xnu);
            printf("TEMPERATURE:  %.3lf\n", E[0]);
            printf("DIAMETER:     %.3lf\n", sigmaVec[0]);
            printf("------------------------------\n");
            break;

        case 4: 
            // DOUBLE YUKAWA (Attractive + Repulsive)
            E[0] = especie1.temperature2/especie1.temperature;  // Attractive strength
            E[2] = especie2.temperature2/especie2.temperature;
            E[1] = sqrt(E[0] * E[2]);

            E2[0] = 1.0/especie1.temperature;  // Repulsive strength
            E2[2] = 1.0/especie2.temperature;
            E2[1] = sqrt(E2[0] * E2[2]);

            z[0] = especie1.lambda;
            z[2] = especie2.lambda;
            z[1] = sqrt(z[0] * z[2]);

            z2[0] = especie1.lambda2;
            z2[2] = especie2.lambda2;
            z2[1] = sqrt(z2[0] * z2[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaVec[k]) {
                        Ua[i*ncols + k] = 0.0;
                        Ur[i*ncols + k] = 0.0;
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        Ua[i*ncols + k] = - E[k] * exp(- z[k] * (r[i] - 1.0)) / r[i];
                        Ur[i*ncols + k] = E2[k] * exp(- z2[k] * (r[i] - 1.0)) / r[i];
                        U[i*ncols + k] = Ua[i*ncols + k] + Ur[i*ncols + k];
                        Up[i*ncols + k] = (1.0 + z[k]*r[i]) * Ua[i*ncols + k] + \
                                          (1.0 + z2[k]*r[i]) * Ur[i*ncols + k];
                    }
                }
            }

            printf("POTENTIAL: DOUBLE YUKAWA (ATRACTIVE + REPULSIVE)\n\n");
            printf("TEMPERATURE  (atr, rep):  %1.9e   %1.9e\n", 1/E[0], 1/E2[0]);
            printf("RATIO  (atr perturbation):  %.3lf\n", E[0]/E2[0]);
            printf("z  (atr, rep):       %.3lf   %.3lf\n", z[0], z2[0]);
            printf("DIAMETER:                 %.3lf\n", sigmaVec[0]);
            printf("------------------------------\n");
            break;

        case 5: 
            // ATTRACTIVE YUKAWA
            E[0] = 1.0 / especie1.temperature;
            E[2] = 1.0 / especie2.temperature;
            E[1] = sqrt(E[0] * E[2]);

            z[0] = especie1.lambda;
            z[2] = especie2.lambda;
            z[1] = sqrt(z[0] * z[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaVec[k]) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        U[i*ncols + k] = - E[k] * exp(- z[k] * (r[i] - 1.0)) / r[i];
                        Up[i*ncols + k] = (1.0 + z[k]*r[i]) * U[i*ncols + k];
                    }
                }
            }

            printf("POTENTIAL: ATRACTIVE YUKAWA\n\n");
            printf("TEMPERATURE:  %.3lf\n", E[0]);
            printf("LAMBDA:       %.3lf\n", z[0]);
            printf("DIAMETER:     %.3lf\n", sigmaVec[0]);
            printf("------------------------------\n");
            break;

        case 6: 
            // REPULSIVE YUKAWA
            E[0] = 1.0 / especie1.temperature;
            E[2] = 1.0 / especie2.temperature;
            E[1] = sqrt(E[0] * E[2]);

            z[0] = especie1.lambda;
            z[2] = especie2.lambda;
            z[1] = sqrt(z[0] * z[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaVec[k]) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        U[i*ncols + k] =  E[k] * exp( -z[k] * (r[i] - 1.0)) / r[i];
                        Up[i*ncols + k] = (1.0 + z[k]*r[i]) * U[i*ncols + k];
                    }
                }
            }
          
            printf("POTENTIAL: REPULSIVE YUKAWA\n\n");
            printf("TEMPERATURE:  %1.9e\n", 1/E[0]);
            printf("LAMBDA:       %1.9e\n", z[0]);
            printf("DIAMETER:     %1.9e\n", sigmaVec[0]);
            printf("------------------------------\n");
            break;

        case 7: 
            // HARD SPHERE
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    U[i*ncols + k] = 0.0;
                    Up[i*ncols + k] = 0.0;
                }
            }

            printf("POTENTIAL: HARD SPHERE\n\n");
            printf("------------------------------\n");
            break;

        case 8: 
            // SHOULDER FUNCTION (Step Potential)
            E[0] = 1.0 / especie1.temperature;
            E[2] = 1.0 / especie2.temperature;
            E[1] = sqrt(E[0] * E[2]);

            z[0] = especie1.lambda;
            z[2] = especie2.lambda;
            z[1] = sqrt(z[0] * z[2]);

            E2[0] = especie1.temperature2;
            E2[2] = especie2.temperature2;
            E2[1] = sqrt(E2[0] * E2[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaVec[k] || (r[i] > E2[k])) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        U[i*ncols + k] =  E[k] * z[k];
                        Up[i*ncols + k] = 0.0;
                    }
                }
            }

            printf("POTENTIAL: SHOULDER FUNCTION\n\n");
            printf("TEMPERATURE:  %.3lf\n", E[0]);
            printf("HEIGHT:       %.3lf\n", z[0]);
            printf("WIDTH:        %.3lf\n", E2[0]);
            printf("DIAMETER:     %.3lf\n", sigmaVec[0]);
            printf("------------------------------\n");
            break;

        case 9: 
            // DOWN-HILL FUNCTION
            E[0] = 1.0 / especie1.temperature;
            E[2] = 1.0 / especie2.temperature;
            E[1] = sqrt(E[0] * E[2]);

            E2[0] = especie1.temperature2;
            E2[2] = especie2.temperature2;
            E2[1] = sqrt(E2[0] * E2[2]);

            z[0] = especie1.lambda / (E2[0] - sigmaVec[0]);
            z[2] = especie2.lambda / (E2[2] - sigmaVec[2]);
            z[1] = sqrt(z[0] * z[2]);
            
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaVec[k] || (r[i] > E2[k])) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        U[i*ncols + k] =  E[k] * (-z[k] * (r[i] - E2[k]));
                        Up[i*ncols + k] = E[k] * z[k] * r[i];
                    }
                }
            }

            printf("POTENTIAL: DOWN-HILL FUNCTION\n\n");
            printf("TEMPERATURE:  %.3lf\n", E[0]);
            printf("HEIGHT:       %.3lf\n", especie1.lambda);
            printf("WIDTH:        %.3lf\n", E2[0]);
            printf("DIAMETER:     %.3lf\n", sigmaVec[0]);
            printf("------------------------------\n");
            break;

        case 10: 
            // GAUSSIAN CORE MODEL
            E[0] = 1.0 / especie1.temperature;
            E[2] = 1.0 / especie2.temperature;
            E[1] = sqrt(E[0] * E[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    U[i*ncols + k] = E[k] * exp(- pow(r[i] / sigmaVec[k], 2.0));
                    Up[i*ncols + k] = 2.0 * pow(r[i] / sigmaVec[k], 2.0) * U[i*ncols + k];
                }
            }

            printf("POTENTIAL: GAUSSIAN CORE MODEL\n\n");
            printf("TEMPERATURE:  %.3lf\n", E[0]);
            printf("DIAMETER:     %.3lf\n", sigmaVec[0]);
            printf("------------------------------\n");
            break;

        case 11: 
            // RAMP
            E[0] = 1.0 / especie1.temperature;
            E[2] = 1.0 / especie2.temperature;
            E[1] = sqrt(E[0] * E[2]);
            double lamb = 1.56;
            for (int k = 0; k < ncols; k++) {
                for (int i = 0; i < nrows; i++) {
                    if (r[i] < sigmaVec[k] || r[i] > lamb * sigmaVec[k]) {
                        U[i*ncols + k]  = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        U[i*ncols + k]  = E[k] * (lamb - r[i]) / (lamb - 1.0);
                        Up[i*ncols + k] = E[k] * r[i] / (lamb - 1.0);
                    }
                }
            }

            printf("POTENTIAL: STEP FUNCTION\n\n");
            printf("TEMPERATURE:  %1.9e\n", 1/E[0]);
            printf("DIAMETER:     %1.9e\n", sigmaVec[0]);
            printf("------------------------------\n");
            break;

        case 12: 
            // STEP FUNCTION (Soft Core)
            E[0] = 1.0 / especie1.temperature;
            E[2] = 1.0 / especie2.temperature;
            E[1] = sqrt(E[0] * E[2]);

            z[0] = especie1.lambda;
            z[2] = especie2.lambda;
            z[1] = sqrt(z[0] * z[2]);

            E2[0] = especie1.temperature2;
            E2[2] = especie2.temperature2;
            E2[1] = sqrt(E2[0] * E2[2]);
            printf("%1.9e\t%1.9e\t%1.9e\n",z[0],E[0],sigmaVec[0]);
            
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] > sigmaVec[k]) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        U[i*ncols + k] =  E[k] * pow(1-r[i]/sigmaVec[0],z[k]); 
                        Up[i*ncols + k] = 0.0;
                    }
                }
            }
            break;

        case 13: 
            // HERTZIAN POTENTIAL (n=2.5)
            // U = E * (1 - r/σ)^n, for r < σ
            double n_exponent = 2.5;

            //E[0] = 1.0; 
            //E[2] = 1.0;
            //E[1] = sqrt(E[0] * E[2]);
            E[0] = 1.0 / especie1.temperature;
            E[2] = 1.0 / especie2.temperature;
            E[1] = sqrt(E[0] * E[2]);
            
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    
                    if (r[i] < sigmaVec[k]) {
                        double factor = 1.0 - (r[i] / sigmaVec[k]);
                        U[i*ncols + k] = E[k] * pow(factor, n_exponent);
                        Up[i*ncols + k] = (n_exponent * E[k] / sigmaVec[k]) * r[i] * pow(factor, n_exponent - 1.0);
                    } else {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    }
                }
            }

            printf("POTENTIAL: HERTZIAN POTENTIAL (n=2.5)\n\n");
            printf("ENERGY SCALE (epsilon/kT):  %.3lf\n", E[0]);
            printf("DIAMETER:                   %.3lf\n", sigmaVec[0]);
            printf("------------------------------\n");
            break;

        case 16:
            // SOFT SHOULDER POTENTIAL
            // U = 0.5 * E * (1 - tanh(alpha * (r - lambda)))
            E[0] = 1.0 / especie1.temperature;
            E[2] = 1.0 / especie2.temperature;
            E[1] = sqrt(E[0] * E[2]);

            z[0]  = especie1.lambda;     // lambda (reach)
            z[2]  = especie2.lambda;
            z[1]  = sqrt(z[0] * z[2]);

            z2[0] = especie1.lambda2;    // alpha (smoothness)
            z2[2] = especie2.lambda2;
            z2[1] = sqrt(z2[0] * z2[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaVec[k]) {
                        U[i*ncols + k]  = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        double arg = z2[k] * (r[i] - z[k]);
                        double t_arg = tanh(arg);
                        U[i*ncols + k]  = 0.5 * E[k] * (1.0 - t_arg);
                        Up[i*ncols + k] = 0.5 * E[k] * z2[k] * r[i] * (1.0 - t_arg * t_arg);
                    }
                }
            }

            printf("POTENTIAL: SOFT SHOULDER POTENTIAL\n\n");
            printf("TEMPERATURE:  %.3lf\n", E[0]);
            printf("LAMBDA (REACH): %.3lf\n", z[0]);
            printf("ALPHA (SMOOTHNESS): %.3lf\n", z2[0]);
            printf("DIAMETER:     %.3lf\n", sigmaVec[0]);
            printf("------------------------------\n");
            break;
    }

    free(Ua);
    free(Ur);
    free(E);
    free(E2);
    free(z);
    free(z2);
}

/**
 * @brief Solves the Ornstein-Zernike equation using Ng's method.
 *
 * This function orchestrates the iterative solution process, handling the
 * charging parameter (kj) for gradual solution and switching between closures.
 *
 * @param Sk Output Structure Factor array.
 * @param Gr Output Radial Distribution Function array.
 * @param potentialID ID of the potential.
 * @param closureID ID of the closure relation (1=PY, 2=HNC, 3=RY).
 * @param alpha Parameter for RY closure.
 * @param EZ Convergence criterion.
 * @param rmax Maximum radial distance.
 * @param nrho Number of density steps.
 * @param folderName Output folder name.
 * @param printFlag Flag to control printing.
 */
void OZ2(double *Sk, double *Gr, int potentialID, int closureID, double alpha, double EZ, \
         double rmax, int nrho, char folderName[20], int *printFlag) {

    int i, k;
    int kj, IRY;
    double T, TFlag, pv, pv0, pv1, pv2;
    double chic, chic0, chic1, chic2;
    double ener, ener0, ener1, ener2;
    double rhoa, dT, drho, ddrho, dalpha, PexV;
    double *cFuncMatrix, *gammaInput1, *gammaInput2, *gammaOutput;
    
    // Allocate memory for solver matrices
    cFuncMatrix = malloc(nrows*ncols * sizeof(double));
    gammaInput1 = malloc(nrows*ncols * sizeof(double));
    gammaInput2 = malloc(nrows*ncols * sizeof(double));
    gammaOutput = malloc(nrows*ncols * sizeof(double));

    if (cFuncMatrix == NULL || gammaInput1 == NULL || gammaInput2 == NULL || gammaOutput == NULL) {
        printf("Memory allocation failed in OZ2.\n");
        return;
    }

    TFlag = 0.0;

    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            gammaInput1[i*ncols + k] = 0.0;
        }
    }

    // Initialize density ramp
    rhoa = rho;
    dT = 1.0 / ((double) nrho);
    drho = rho / ((double) nrho);

    kj = 1;
    rho = kj * drho;
    T = dT * kj;

    // Initial guess
    Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, TFlag, alpha, EZ, rmax, nrho, printFlag);

    // Ramp up density
    while (kj <= 1) {
        for (k = 0; k < ncols; k++) {
            for (i = 0; i < nrows; i++) {
                gammaInput1[i*ncols + k] = gammaOutput[i*ncols + k];
            }
        }
        kj++;
        rho = kj * drho;
        T = dT * kj;
        Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, TFlag, alpha, EZ, rmax, nrho, printFlag);
    }

    Extrap(gammaInput1, gammaOutput, rho, drho);

    while(true) {
        for (k = 0; k < ncols; k++) {
            for (i = 0; i < nrows; i++) {
                gammaInput2[i*ncols + k] = gammaOutput[i*ncols + k];
            }
        }

        kj++;
        rho = kj * drho;
        T = dT * kj;

        Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, TFlag, alpha, EZ, rmax, nrho, printFlag);

        if (kj == nrho) {
            break;
        } else {
            Extrap(gammaInput2, gammaOutput, rho, drho);
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    gammaInput1[i*ncols + k] = gammaInput2[i*ncols + k];
                }
            }
        }
    }

    // Final solution step
    T = 1.0;
    rho = rhoa;

    Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, TFlag, alpha, EZ, rmax, nrho, printFlag);
    
    switch (closureID){
        case 1:
            printf("\n==============\n");
            printf("  Salida PY:\n");
            printf("==============\n\n");
            break;
        case 2:
            printf("\n==============\n");
            printf(" Salida HNC:\n");
            printf("==============\n\n");
            break;
        case 3:
            // Rogers-Young specific logic to determine alpha
            printf("\n===========================\n");
            printf("  CALCULANDO VALOR ALPHA:\n");
            printf("===========================\n\n");

            ddrho = rhoa / 100.0;
            dalpha = -alpha / 50.0;

            do{
                T = 1.0;
                rho = rhoa;

                Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, TFlag, alpha, EZ, rmax, nrho, printFlag);
                
                for (k = 0; k < ncols; k++) {
                    for (i = 0; i < nrows; i++) {
                        gammaInput1[i*ncols + k] = gammaOutput[i*ncols + k];
                    }
                }

                Termo(gammaOutput, cFuncMatrix, &pv0, &chic0, &ener0);
                chic = chic0;

                rho -= ddrho;
                Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, TFlag, alpha, EZ, rmax, nrho, printFlag);
                Termo(gammaOutput, cFuncMatrix, &pv1, &chic1, &ener1);

                rho += 2.0*ddrho;
                Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, TFlag, alpha, EZ, rmax, nrho, printFlag);
                Termo(gammaOutput, cFuncMatrix, &pv2, &chic2, &ener2);

                RY(pv1, pv2, chic, ddrho, &alpha, dalpha, &IRY);

            } while(IRY == 1);
            
            printf("\n==============\n");
            printf("  Salida RY:\n");
            printf("==============\n\n");
            break;
    }

    // Final calculation and output
    T = 1.0;
    TFlag = 1.0;
    rho = rhoa;

    Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, TFlag, alpha, EZ, rmax, nrho, printFlag);

    Escribe(gammaOutput, cFuncMatrix, Sk, Gr, potentialID, closureID, folderName);

    Termo(gammaOutput, cFuncMatrix, &pv, &chic, &ener);
    PexV = pv/rho - 1.0;

    free(cFuncMatrix);
    free(gammaInput1);
    free(gammaInput2);
    free(gammaOutput);
}

/**
 * @brief Calculates thermodynamic properties (Pressure, Compressibility, Energy).
 *
 * @param gamma Indirect correlation function.
 * @param cFuncMatrix Direct correlation function.
 * @param pv1 Pointer to store pressure (virial).
 * @param chic Pointer to store compressibility.
 * @param ener Pointer to store energy.
 */
void Termo(double *gamma, double *cFuncMatrix, double *pv1, double *chic, double *ener) {
    
    int i, k;
    double ru1, ru2, ru3;
    double *r1;
    double *gMatrix;

    r1 = malloc(nrows * sizeof(double));
    gMatrix = malloc(nrows*ncols * sizeof(double));

    if (r1 == NULL || gMatrix == NULL) {
        printf("Memory allocation failed in Termo.\n");
        return;
    }

    for (i = 0; i < nrows; i++) {
        r1[i] = x[0]*x[0]*cFuncMatrix[i*ncols + 0] + 2.0*x[0]*x[1]*cFuncMatrix[i*ncols + 1] + x[1]*x[1]*cFuncMatrix[i*ncols + 2];
        r1[i] = r1[i] * r[i]*r[i];
    }

    *chic = 0.0;
    
    for (i = 1; i < nrows - 1; i++) {
        *chic += r1[i];
    }
    *chic = dr * (*chic + (r1[0] + r1[nrows-1]) / 2.0);
    *chic = 1.0 - 4.0 * M_PI * rho * (*chic);

    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            gMatrix[i + k*nrows] = gamma[i*ncols + k] + cFuncMatrix[i*ncols + k] + 1.0;
        }
    }

    for (i = 0; i < nrows; i++) {
        ru1 = x[0]*x[0] * gMatrix[i + 0*nrows] * Up[i*ncols + 0];
        ru2 = 2.0 * x[0]*x[1] * gMatrix[i + 1*nrows] * Up[i*ncols + 1];
        ru3 = x[1]*x[1] * gMatrix[i + 2*nrows] * Up[i*ncols + 2];
        r1[i] = (ru1 + ru2 + ru3) * r[i]*r[i];
    }

    *pv1 = 0.0;

    for (i = 1; i < (nrows-1); i++) {
        *pv1 += r1[i];
    }
    
    *pv1 = dr * ((*pv1) + (r1[0] + r1[nrows-1])/2.0);
    *pv1 = rho * (1.0 + 2.0*M_PI * rho*(*pv1)/3.0);

    for (i = 0; i < nrows; i++) {
        r1[i] = x[0]*x[0] * gMatrix[i + 0*nrows] * U[i*ncols + 0] + 2.0*x[0]*x[1] * gMatrix[i + 1*nrows]*U[i*ncols + 1];
        r1[i] = (r1[i] + x[1]*x[1] * gMatrix[i + 2*nrows] * U[i*ncols + 2]) * r[i]*r[i];
    }

    *ener = 0.0;

    for (i = 1; i < (nrows-1); i++) {
        *ener += r1[i];
    }

    *ener = dr * ((*ener) + (r1[0] + r1[nrows-1]) / 2.0);
    *ener *= 2.0 * M_PI * rho;

    free(r1);
    free(gMatrix);
}

/**
 * @brief Iteratively adjusts alpha for the Rogers-Young closure.
 *
 * Checks for thermodynamic consistency between virial and compressibility equations of state.
 *
 * @param pv1 Pressure at rho - drho.
 * @param pv2 Pressure at rho + 2*drho.
 * @param chic Compressibility.
 * @param ddrho Density step.
 * @param alpha Pointer to alpha parameter.
 * @param dalpha Alpha step.
 * @param IRY Pointer to flag (1 if iteration needed, 0 if converged).
 */
void RY(double pv1, double pv2, double chic, double ddrho, double *alpha, double dalpha, int *IRY) {

    static double dif[2] = {0.0};
    static int ix = 1;
    double chiv, prod, A, B;

    *IRY = 0;
    chiv = (pv2 - pv1) / (2.0 * ddrho);
    dif[ix-1] = chic - chiv;

    prod = 0.0;

    if (ix < 2) {
        ix++;
        *alpha += dalpha;
        *IRY = 1;
        printf("   CHIC = %.17g   CHIV = %.17g\n", chic, chiv);
        printf("   DIFF = %.17g   \n", dif[ix - 1]);
        printf("   ALPHA = %.17g  <------------------ \n\n", *alpha);
        return;
    }

    prod = dif[0] * dif[1];

    if (prod < 0.0) {
        *alpha -= dalpha;
        A = (dif[1] - dif[0]) / dalpha;
        B = dif[0] - A * (*alpha);
        *alpha = -B / A;
        *IRY = 0;
        printf("   CHIC = %.17g   CHIV = %.17g\n", chic, chiv);
        printf("   DIFF = %.17g   \n", dif[ix - 1]);
        printf("   ALPHA = %.17g  <------------------ \n\n", *alpha);
        return;
    } else {
        dif[0] = dif[1];
        *alpha += dalpha;
        *IRY = 1;
        printf("   CHIC = %.17g   CHIV = %.17g\n", chic, chiv);
        printf("   DIFF = %.17g   \n", dif[ix - 1]);
        printf("   ALPHA = %.17g  <------------------ \n\n", *alpha);
        return;
    }
}

/**
 * @brief Writes the results (Structure Factor and Radial Distribution Function) to files.
 *
 * @param gamma Indirect correlation function.
 * @param cFuncMatrix Direct correlation function.
 * @param Sk Output array for Structure Factor.
 * @param Gr Output array for Radial Distribution Function.
 * @param potentialID Potential ID.
 * @param closureID Closure ID.
 * @param folderName Output folder name.
 */
void Escribe(double *gamma, double *cFuncMatrix, double *Sk, double *Gr, int potentialID, int closureID, char folderName[20]) {

    int i, k;
    double dk, qmax, rk_max, sqmax, delta;
    double *rk, *c1, *gh, *gh2, *Ck, *S;

    rk  = malloc(nrows * sizeof(double));
    c1  = malloc(nrows*ncols * sizeof(double));
    gh  = malloc(nrows*ncols * sizeof(double));
    gh2 = malloc(nrows*ncols * sizeof(double));
    Ck  = malloc(nrows*ncols * sizeof(double));
    S   = malloc(nrows*ncols * sizeof(double));

    if (rk == NULL || c1 == NULL || gh == NULL || gh2 == NULL || Ck == NULL || S == NULL) {
        printf("Memory allocation failed in Escribe.\n");
        return;
    }

    for (i = 0; i < nrows; i++) {
        for (k = 0; k < ncols; k++) {
            gh[i*ncols + k] = gamma[i*ncols + k] + cFuncMatrix[i*ncols + k];
        }
        Gr[i*2 + 0] = r[i];
        Gr[i*2 + 1] = gh[i*ncols + 0] + 1.0;
    }

    qmax = q[nrows - 1];
    rk_max = qmax / 2.0;
    dk = rk_max / (1.0 * nrows);

    for (i = 0; i < nrows; i++) {
        rk[i] = 1.0E-5 + (i+0) * dk;
    }

    FT(gh, gh2, rk, dr);
    FT(cFuncMatrix, c1, rk, dr);

    for (i = 0; i < nrows; i++) {
        for (k = 0; k < ncols; k++) {
            Ck[i*ncols + k] = c1[i*ncols + k];
        }
    }

    for (i = 0; i < nrows; i++) {
        delta = (1.0 - rho*x[0] * Ck[i*ncols + 0]) * (1.0 - rho*x[1] * Ck[i*ncols + 2]);
        delta -= pow(rho, 2.0) * x[0]*x[1] * pow(Ck[i*ncols + 1], 2.0);
        
        S[i*ncols + 0] = (1.0 - rho*x[1] * Ck[i*ncols + 2]) * Ck[i*ncols + 0];
        S[i*ncols + 0] = (S[i*ncols + 0] + rho*x[1] * pow(Ck[i*ncols + 1], 2.0)) / delta;
        S[i*ncols + 1] = Ck[i*ncols + 1] / delta;
        S[i*ncols + 2] = (1.0 - rho*x[0] * Ck[i*ncols + 0]) * Ck[i*ncols + 2];
        S[i*ncols + 2] = (S[i*ncols + 2] + rho*x[0] * pow(Ck[i*ncols + 1], 2.0)) / delta;
        S[i*ncols + 0] = x[0] + rho*pow(x[0], 2.0) * S[i*ncols + 0];
        S[i*ncols + 1] = rho * x[0]*x[1] * S[i*ncols + 1];
        S[i*ncols + 2] = x[1] + rho*pow(x[1], 2.0) * S[i*ncols + 2];
    }

    sqmax = 0.0;

    for (i = 0; i < nrows; i++) {
        Sk[i*2 + 0] = rk[i];
        Sk[i*2 + 1] = S[i*ncols + 0]/x[0];
        if (S[i*ncols + 0] > sqmax) {
            sqmax = S[i*ncols + 0];
        }
    }

    free(rk);
    free(c1);
    free(gh);
    free(gh2);
    free(Ck);
    free(S);
}

/**
 * @brief Ng's method for accelerating convergence of the iterative solution.
 *
 * @param kj Current iteration step (related to density ramp).
 * @param gammaInput Input gamma function.
 * @param gammaOutput Output gamma function.
 * @param potentialID Potential ID.
 * @param closureID Closure ID.
 * @param cFuncMatrix Direct correlation function.
 * @param T Temperature parameter.
 * @param TFlag Temperature flag.
 * @param alpha Alpha parameter.
 * @param EZ Convergence criterion.
 * @param rmax Maximum radial distance.
 * @param nrho Number of density steps.
 * @param printFlag Print flag.
 */
void Ng(int kj, double *gammaInput, double *gammaOutput, int potentialID, int closureID, double *cFuncMatrix, \
        double T, double TFlag, double alpha, double EZ, double rmax, int nrho, int *printFlag) {

    int i, k, flag;
    double ETA, V, condition1;
    double *f, *g1, *g2, *g3;
    double *d1, *d2, *d3, *d01, *d02;
    double *d01d01, *d01d02, *d02d02, *d3d01, *d3d02;
    double *const1, *const2;

    f   = malloc(nrows*ncols * sizeof(double));
    g1  = malloc(nrows*ncols * sizeof(double));
    g2  = malloc(nrows*ncols * sizeof(double));
    g3  = malloc(nrows*ncols * sizeof(double));
    d1  = malloc(nrows*ncols * sizeof(double));
    d2  = malloc(nrows*ncols * sizeof(double));
    d3  = malloc(nrows*ncols * sizeof(double));
    d01 = malloc(nrows*ncols * sizeof(double));
    d02 = malloc(nrows*ncols * sizeof(double));
    d01d01  = malloc(ncols * sizeof(double));
    d01d02  = malloc(ncols * sizeof(double));
    d02d02  = malloc(ncols * sizeof(double));
    d3d01   = malloc(ncols * sizeof(double));
    d3d02   = malloc(ncols * sizeof(double));
    const1  = malloc(ncols * sizeof(double));
    const2  = malloc(ncols * sizeof(double));

    if (f == NULL || g1 == NULL || g2 == NULL || g3 == NULL || d1 == NULL || d2 == NULL || d3 == NULL ||
        d01 == NULL || d02 == NULL || d01d01 == NULL || d01d02 == NULL || d02d02 == NULL ||
        d3d01 == NULL || d3d02 == NULL || const1 == NULL || const2 == NULL ) {
        printf("Memory allocation failed in Ng.\n");
        return;
    }

    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            f[i*ncols + k] = gammaInput[i*ncols + k];
        }
    }

    ONg(f, g1, potentialID, closureID, cFuncMatrix, T, TFlag, alpha, rmax);
    ONg(g1, g2, potentialID, closureID, cFuncMatrix, T, TFlag, alpha, rmax);
    ONg(g2, g3, potentialID, closureID, cFuncMatrix, T, TFlag, alpha, rmax);

    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            d1[i*ncols + k] = (g1[i*ncols + k] - f[i*ncols + k]);
            d2[i*ncols + k] = (g2[i*ncols + k] - g1[i*ncols + k]);
        }
    }

    while (kj >= 2) {
        
        for (k = 0; k < ncols; k++) {
            for (i = 0; i < nrows; i++) {
                d3[i*ncols + k] = (g3[i*ncols + k] - g2[i*ncols + k]);
            }
        }

        for (k = 0; k < ncols; k++) {
            for (i = 0; i < nrows; i++) {
                d01[i*ncols + k] = (d3[i*ncols + k] - d2[i*ncols + k]);
                d02[i*ncols + k] = (d3[i*ncols + k] - d1[i*ncols + k]);
            }
        }

        pp(dr, d01, d01, d01d01);
        pp(dr, d01, d02, d01d02);
        pp(dr, d3, d01, d3d01);
        pp(dr, d02, d02, d02d02);
        pp(dr, d3, d02, d3d02);

        V = (double) 1.0E-50;

        for (k = 0; k < ncols; k++) {

            condition1 = d02d02[k] - (d01d02[k]*d01d02[k]) / d01d01[k];

            if ( ((d01d01[k] > V) && (d01d02[k] > V)) && (condition1 > V) ) {
                    
                const2[k] = d3d02[k] - d01d02[k]*d3d01[k] / d01d01[k];
                const2[k] = const2[k] / (d02d02[k] - d01d02[k]*d01d02[k] / d01d01[k]);
                const1[k] = (d3d02[k] - d02d02[k] * const2[k]) / d01d02[k];
                
                for (i = 0; i < nrows; i++) {
                    f[i*ncols + k] = (1.0 - const1[k] - const2[k]) * g3[i*ncols + k];
                    f[i*ncols + k] = f[i*ncols + k] + const1[k]*g2[i*ncols + k] + const2[k]*g1[i*ncols + k];
                }

                flag = 0;

            } else {          
                flag = 1;
                break;
            }
        }

        if (flag == 0){
            ONg(f, g3, potentialID, closureID, cFuncMatrix, T, TFlag, alpha, rmax);
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    d3[i*ncols + k] = g3[i*ncols + k] - f[i*ncols + k];
                }
            }
        }

        Pres(d3, dr, &ETA);

        if (ETA <= EZ) {
            break;
        }

        for (k = 0; k < ncols; k++) {
            for (i = 0; i < nrows; i++) {
                g1[i*ncols + k] = g2[i*ncols + k];
                d1[i*ncols + k] = d2[i*ncols + k];
                g2[i*ncols + k] = g3[i*ncols + k];
                d2[i*ncols + k] = d3[i*ncols + k];
            }
        }

        ONg(g2, g3, potentialID, closureID, cFuncMatrix, T, TFlag, alpha, rmax);
    }

    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            gammaOutput[i*ncols + k] = g3[i*ncols + k];
        }
    }

    if ((*printFlag) == 0){
        printf("%03d/%03d", kj, nrho);
        fflush(stdout);
        printf("\b\b\b\b\b\b\b");
    }

    if (kj >= nrho){
        *printFlag = 1;
    }
    
    free(f);
    free(g1);
    free(g2);
    free(g3);
    free(d1);
    free(d2);
    free(d3);
    free(d01);
    free(d02);
    free(d01d01);
    free(d01d02);
    free(d02d02);
    free(d3d01);
    free(d3d02);
    free(const1);
    free(const2);
}

/**
 * @brief Applies the closure relation (PY, HNC, RY) to calculate the direct correlation function.
 *
 * @param gamma Indirect correlation function.
 * @param potentialID Potential ID.
 * @param closureID Closure ID.
 * @param cFuncMatrix Output Direct correlation function.
 * @param T Temperature parameter.
 * @param alpha Alpha parameter (for RY).
 */
void closrel(double *gamma, int potentialID, int closureID, double *cFuncMatrix, double T, double alpha) {
    
    int i, k;
    double sigmaAux[ncols];
    double arg, F;

    for (k = 0; k < ncols; k++) {
        if (potentialID == 1 || potentialID == 2 || potentialID == 3) {
            sigmaAux[k] = (sigmaVec[k] / 2.0);
        } else if (potentialID == 10) {
            sigmaAux[k] = 0.0;
        } else {
            sigmaAux[k] = sigmaVec[k];
        }
    }

    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            cFuncMatrix[i*ncols + k] = gamma[i*ncols + k] + 1.0;
        }
    }

    switch(closureID){
        
        case 1: // PY
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaAux[k]) {
                        cFuncMatrix[i*ncols + k] = -cFuncMatrix[i*ncols + k];
                    } else {
                        arg = U[i*ncols + k] * T;
                        if (arg > 70.0) {
                            cFuncMatrix[i*ncols + k] = -cFuncMatrix[i*ncols + k];
                        } else {
                            cFuncMatrix[i*ncols + k] = -(exp(-arg) - 1.0) * cFuncMatrix[i*ncols + k];
                        }
                    }
                }
            }
            return;
        
        case 2: // HNC
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaAux[k]) {
                        cFuncMatrix[i*ncols + k] = -cFuncMatrix[i*ncols + k];
                    } else {
                        arg = U[i*ncols + k] * T - gamma[i*ncols + k];
                        if (arg > 70.0) {
                            cFuncMatrix[i*ncols + k] = -cFuncMatrix[i*ncols + k];
                        } else {
                            cFuncMatrix[i*ncols + k] = exp(-arg) - cFuncMatrix[i*ncols + k];
                        }
                    }
                }
            }
            return;
        
        case 3: // RY
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaAux[k]) {
                        cFuncMatrix[i*ncols + k] = -cFuncMatrix[i*ncols + k];
                    } else {
                        F = 1.0 - exp(-alpha * r[i]);
                        arg = (exp(gamma[i*ncols + k] * F) - 1.0) / F;
                        cFuncMatrix[i*ncols + k] = exp(-U[i*ncols + k] * T) * (1.0 + arg) - cFuncMatrix[i*ncols + k];
                    }
                }
            }
            return;
    }
}

/**
 * @brief Appends the closure name to a string.
 *
 * @param inputString String to append to.
 * @param closureID Closure ID.
 */
void appendclosureID(char *inputString, int closureID) {
    switch (closureID){
        case 1: strcat(inputString, "_PY"); break;
        case 2: strcat(inputString, "_HNC"); break;
        case 3: strcat(inputString, "_RY"); break;
    }
}

/**
 * @brief Appends the potential name to a string.
 *
 * @param inputString String to append to.
 * @param potentialID Potential ID.
 */
void appendPotentialID(char *inputString, int potentialID) {
    switch (potentialID){
        case 1: strcat(inputString, "_IPL"); break;
        case 2: strcat(inputString, "_LJT"); break;
        case 3: strcat(inputString, "_LJT2"); break;
        case 4: strcat(inputString, "_DBLEYUK"); break;
        case 5: strcat(inputString, "_ATRYUK"); break;
        case 6: strcat(inputString, "_REPYUK"); break;
        case 7: strcat(inputString, "_HS"); break;
        case 8: strcat(inputString, "_STEPFUNC"); break;
        case 9: strcat(inputString, "_DOWNHILL"); break;
        case 16: strcat(inputString, "_SOFTSHOULDER"); break;
    }
}

/**
 * @brief Prints the name of the potential.
 *
 * @param potentialID Potential ID.
 * @param xnu Potential parameter.
 */
void PotentialName(int potentialID, double xnu) {
    switch (potentialID){
        case 1: printf("POTENTIAL: INVERSE POWER LAW"); break;
        case 2: printf("POTENTIAL: TRUNCATED LENNARD-JONES 6-12"); break;
        case 3: printf("POTENTIAL: TRUNCATED LENNARD-JONES %.1lf-%.1lf", xnu, 2.0*xnu); break;
        case 4: printf("POTENTIAL: DOUBLE YUKAWA (ATRACTIVE + REPULSIVE)"); break;
        case 5: printf("POTENTIAL: ATRACTIVE YUKAWA"); break;
        case 6: printf("POTENTIAL: REPULSIVE YUKAWA"); break;
        case 7: printf("POTENTIAL: HARD SPHERE"); break;
        case 8: printf("POTENTIAL: SHOULDER FUNCTION"); break;
        case 9: printf("POTENTIAL: DOWN-HILL FUNCTION"); break;
        case 16: printf("POTENTIAL: SOFT SHOULDER POTENTIAL"); break;
    }
}
