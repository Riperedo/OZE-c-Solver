#ifndef FACDES2Y_FUNCS_DOT_H    /* This is an "include guard" */
#define FACDES2Y_FUNCS_DOT_H    /* prevents the file from being included twice. */

#include "math_aux.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

// We define some global parameters
extern int nrows, ncols;
extern double rho, dr;
extern double *r, *q, x[2];
extern double *U, *Up, *sigmaVec;

typedef struct species{

    double diameter;
    double temperature;
    double lambda;
    double temperature2; //Used if we are considering a double Yukawa potential
    double lambda2;      //Used if we are considering a double Yukawa potential

}species;

void input(double fv, double xnu, species especie1, species especie2, double rmax, int potentialID);
void POT(species especie1, species especie2, int potentialID, double xnu);
void OZ2(double *Sk, double *Gr, double *Gamma, int potentialID, int closureID, double alpha, double EZ, \
         double rmax, int nrho, char folderName[20], int *printFlag);
void Termo(double *gamma, double *cFuncMatrix, double *pv1, double *chic, double *ener);
void RY(double pv1, double pv2, double chic, double ddrho, double *alpha, double dalpha, int *IRY);
void Escribe(double *gamma, double *cFuncMatrix, double *Sk, double *Gr, double *Gamma, int potentialID, int closureID, char folderName[20]);
void Ng(int kj, double *gammaInput, double *gammaOutput, int potentialID, int closureID, double *cFuncMatrix, \
        double T, double Tfin, double alpha, double EZ, double rmax, int nrho, int *printFlag); 
void closrel(double *gamma, int potentialID, int closureID, double *cFuncMatrix, double T, double alpha);

void appendclosureID(char *inputString, int closureID);
void appendPotentialID(char *inputString, int potentialID);
void PotentialName(int potentialID, double xnu);
void printLoadingBar(int progress, int total);

void set_init_gamma_file(const char *filename);
const char *get_init_gamma_file(void);
int load_init_gamma(const char *filename, const double *r_grid, double *gammaInput, int nrows, int ncols);

#endif /* FACDES2Y_FUNCS_DOT_H */
