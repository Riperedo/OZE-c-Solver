#ifndef FACDES2Y_DOT_H    /* This is an "include guard" */
#define FACDES2Y_DOT_H    /* prevents the file from being included twice. */

#include "structures.h"
#include "math_aux.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

void interpolationFunc(double *xInput, double *yInput, double *xOutput, double *yOutput, int nrowsInput, int nrowsOutput);

void ck_HNC(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y);
void is_HNC(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y);
void sk_HNC(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y);
void gr_HNC(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *r, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y);

void ck_PY(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y);
void is_PY(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y);
void sk_PY(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y);
void gr_PY(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *r, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y);

void ck_RY(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y);
void is_RY(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y);
void sk_RY(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y);
void gr_RY(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *r, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y);


int facdes2YFunc(const int nodes, int nrho, double rmax, int potentialID, int closureID, double sigma1, double sigma2, \
                 double Temperature, double Temperature2, double lambda_a, double lambda_r, double volumeFactor, \
                 double d, double alpha, double EZ, const int OutputFlag, double *ykVec, double *rkVec);

char *getFolderID();
bool directoryExists(char *path);

#endif /* FACDES2Y_DOT_H */
