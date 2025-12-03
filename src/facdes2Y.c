/**
 * @file facdes2Y.c
 * @brief Interface for solving the Ornstein-Zernike equation using HNC and Rogers-Young closures.
 *
 * This file contains functions to calculate the Direct Correlation Function,
 * Inverse Structure Factor, and Structure Factor for colloidal systems.
 * It uses the Hypernetted Chain (HNC) and Rogers-Young (RY) approximations.
 */

#include "facdes2Y.h"
#include <sys/stat.h> // For mkdir

/**
 * @brief Maximum range for the radial distribution function.
 */
double rmax;

/**
 * @brief Number of rows and columns for data structures.
 */
int nrows, ncols;

/**
 * @brief Number of density points.
 */
int nrho = 100;

/**
 * @brief Diameter scaling factor.
 */
double d = 1.0;

/**
 * @brief Potential parameter xnu.
 */
double xnu = 14.0;

/**
 * @brief Alpha parameter for the closure.
 */
double alpha = 1.0;

/**
 * @brief EZ parameter.
 */
double EZ = 1.0E-4;

/**
 * @brief Diameter of species 1.
 */
double sigma1 = 1.0;

/**
 * @brief Diameter of species 2.
 */
double sigma2 = 1.0;

/**
 * @brief Global density and spatial step parameters.
 */
double rho, dr;

/**
 * @brief Arrays for position, momentum, and mole fractions.
 */
double *r, *q, x[2];

/**
 * @brief Potential energy arrays and sigma vector.
 */
double *U, *Up, *sigmaVec;

// Helper function prototype
static void solve_and_process(double volumeFactor, double Temperature, double Temperature2, 
                            double lambda_a, double lambda_r, const gsl_vector *input_vec, 
                            double *OutputVec, int potentialNumber, int nodesFacdes2Y,
                            int closureID, int outputFlag, const char *filename);

/**
 * @brief Calculates the Direct Correlation Function using the HNC closure.
 */
void ck_HNC(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y){
    solve_and_process(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, k, OutputVec, 
                     potentialNumber, nodesFacdes2Y, 2, 1, "HNC_CdeK.dat");
}

/**
 * @brief Calculates the Inverse Structure Factor using the HNC closure.
 */
void is_HNC(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y){
    solve_and_process(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, k, OutputVec, 
                     potentialNumber, nodesFacdes2Y, 2, 2, "HNC_FT_CdeK.dat");
}

/**
 * @brief Calculates the Structure Factor using the HNC closure.
 */
void sk_HNC(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y){
    solve_and_process(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, k, OutputVec, 
                     potentialNumber, nodesFacdes2Y, 2, 0, "HNC_SdeK.dat");
}

/**
 * @brief Calculates the Direct Correlation Function using the Rogers-Young closure.
 */
void ck_RY(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y){
    solve_and_process(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, k, OutputVec, 
                     potentialNumber, nodesFacdes2Y, 3, 1, "RY_CdeK.dat");
}

/**
 * @brief Calculates the Inverse Structure Factor using the Rogers-Young closure.
 */
void is_RY(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y){
    solve_and_process(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, k, OutputVec, 
                     potentialNumber, nodesFacdes2Y, 3, 2, "RY_FT_CdeK.dat");
}

/**
 * @brief Calculates the Structure Factor using the Rogers-Young closure.
 */
void sk_RY(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *k, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y){
    solve_and_process(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, k, OutputVec, 
                     potentialNumber, nodesFacdes2Y, 3, 0, "RY_SdeK.dat");
}

/**
 * @brief Calculates the Radial Distribution Function using the HNC closure.
 */
void gr_HNC(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *r_vec, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y){
    solve_and_process(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, r_vec, OutputVec, 
                     potentialNumber, nodesFacdes2Y, 2, 3, "HNC_GdeR.dat");
}

/**
 * @brief Calculates the Radial Distribution Function using the Rogers-Young closure.
 */
void gr_RY(double volumeFactor, double Temperature, double Temperature2, double lambda_a, double lambda_r, const gsl_vector *r_vec, \
            double *OutputVec, int potentialNumber, int nodesFacdes2Y){
    solve_and_process(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, r_vec, OutputVec, 
                     potentialNumber, nodesFacdes2Y, 3, 3, "RY_GdeR.dat");
}

/**
 * @brief Generic helper function to solve OZ equation and process results.
 * 
 * This function encapsulates the common logic for all solver variants:
 * 1. Memory allocation
 * 2. Solving the OZ equation via facdes2YFunc
 * 3. Interpolating results to the input grid
 * 4. Writing results to file in the output/ directory
 * 5. Cleaning up memory
 */
static void solve_and_process(double volumeFactor, double Temperature, double Temperature2, 
                            double lambda_a, double lambda_r, const gsl_vector *input_vec, 
                            double *OutputVec, int potentialNumber, int nodesFacdes2Y,
                            int closureID, int outputFlag, const char *filename) {
    
    int potentialID = potentialNumber;
    const int nodesInput = input_vec->size;
    
    double *rkVec, *ykVec, *xOutputVec;
    
    rmax = 160;
    
    // Allocate memory for solver arrays and output
    rkVec = malloc(nodesFacdes2Y * sizeof(double));
    ykVec = malloc(nodesFacdes2Y * sizeof(double));
    xOutputVec = malloc(nodesInput * sizeof(double));
    
    if (rkVec == NULL || ykVec == NULL || xOutputVec == NULL) {
        printf("Memory allocation failed in solve_and_process for %s.\n", filename);
        if (rkVec) free(rkVec);
        if (ykVec) free(ykVec);
        if (xOutputVec) free(xOutputVec);
        return;
    }
    
    for(int i=0; i < nodesInput; i++){
        xOutputVec[i] = (double)(input_vec->data[i]);
        OutputVec[i] = 0.0;
    }
    
    // Solve OZ equation
    facdes2YFunc(nodesFacdes2Y, nrho, rmax, potentialID, closureID, sigma1, sigma2, Temperature, Temperature2, lambda_a, lambda_r, \
                 volumeFactor, d, alpha, EZ, outputFlag, ykVec, rkVec);
    
    // Interpolate results to input grid
    interpolationFunc(rkVec, ykVec, xOutputVec, OutputVec, nodesFacdes2Y, (int)nodesInput);
    
    // Write to file in output/ directory
    char filepath[256];
    snprintf(filepath, sizeof(filepath), "output/%s", filename);
    
    FILE *outputFile = fopen(filepath, "w");
    if (outputFile == NULL) {
        // Try creating the directory if it doesn't exist (though Makefile should handle this)
        // Or fallback to current directory if output/ fails
        fprintf(stderr, "Warning: Could not open %s for writing. Trying local directory.\n", filepath);
        outputFile = fopen(filename, "w");
    }
    
    if (outputFile != NULL) {
        for (int i=0; i < nodesFacdes2Y; i++) {
            fprintf(outputFile, "%.17lf\t%.17lf\n", rkVec[i], ykVec[i]);
            // For S(k) HNC, the original code also printed to stdout. Preserving that behavior if needed, 
            // but usually libraries shouldn't print to stdout. Removed for consistency unless requested.
            // if (outputFlag == 0 && closureID == 2) printf("%.17lf\t%.17lf\n", rkVec[i], ykVec[i]);
        }
        fclose(outputFile);
    } else {
        fprintf(stderr, "Error: Could not open file %s for writing.\n", filename);
    }
    
    free(xOutputVec);
    free(rkVec);
    free(ykVec);
}

/**
 * @brief Main solver function for the Ornstein-Zernike equation.
 *
 * @param nodes Number of spatial nodes.
 * @param nrho Number of density points.
 * @param rmax Maximum range.
 * @param potentialID Interaction potential ID.
 * @param closureID Closure relation ID.
 * @param sigma1 Diameter of species 1.
 * @param sigma2 Diameter of species 2.
 * @param Temperature Temperature of species 1.
 * @param Temperature2 Temperature of species 2.
 * @param lambda_a Attraction range.
 * @param lambda_r Repulsion range.
 * @param volumeFactor Volume fraction.
 * @param d Diameter scaling.
 * @param alpha Closure parameter.
 * @param EZ EZ parameter.
 * @param OutputFlag Flag to determine output type.
 * @param ykVec Output array for function values.
 * @param rkVec Output array for radial/wave vector values.
 * @return 0 on success, 1 on failure.
 */
int facdes2YFunc(const int nodes, int nrho, double rmax, int potentialID, int closureID, double sigma1, double sigma2, \
                 double Temperature, double Temperature2, double lambda_a, double lambda_r, double volumeFactor, \
                 double d, double alpha, double EZ, const int OutputFlag, double *ykVec, double *rkVec) {
    
    nrows = nodes;
    ncols = 3;

    int i;
    int printFlag = 0;
    double *StructFactor, *FT_Cr, *Gr_data;
    bool IsPolidispersed;
    species especie1, especie2;
    
    // Allocate memory for global arrays
    r               = malloc(nrows * sizeof(double));
    q               = malloc(nrows * sizeof(double));
    U               = malloc(nrows*ncols * sizeof(double));
    Up              = malloc(nrows*ncols *sizeof(double));
    sigmaVec        = malloc(ncols * sizeof(double));
    StructFactor    = malloc(nrows*2 * sizeof(double));
    FT_Cr           = malloc(nrows*2 * sizeof(double));
    Gr_data         = malloc(nrows*2 * sizeof(double));

    if (r == NULL || q == NULL || U == NULL || Up == NULL || sigmaVec == NULL || StructFactor == NULL || FT_Cr == NULL || Gr_data == NULL) {
        printf("Memory allocation failed in facdes2YFunc.\n");
        return 1;
    }
    
    x[0] = 1.0;
    x[1] = 1.0 - x[0];

    IsPolidispersed = false;

    especie1.diameter = sigma1;
    especie1.temperature = Temperature;
    especie1.lambda = lambda_a;
    especie1.temperature2 = Temperature2;
    especie1.lambda2 = lambda_r;

    if (IsPolidispersed) {
        especie2.diameter = sigma2;
        especie2.temperature = Temperature;
        especie2.lambda = lambda_a;
        especie2.temperature2 = Temperature2;
        especie2. lambda2 = lambda_r;
    } else {
        especie2.diameter = especie1.diameter;
        especie2.temperature = especie1.temperature;
        especie2.lambda = especie1.lambda;
        especie2.temperature2 = especie1.temperature2;
        especie2. lambda2 = especie1. lambda2;
    }

    char *folderName = getFolderID();

    // Read input data
    input(volumeFactor, xnu, especie1, especie2, rmax, potentialID);

    // Perform calculations
    OZ2(StructFactor, Gr_data, potentialID, closureID, alpha, EZ, rmax, nrho, folderName, &printFlag);

    printf("\n\n");

    switch(OutputFlag){
        
        case 1: // Fourier transform HNC closure
            for (i=0; i<nrows; i++){
                rkVec[i] = FT_Cr[i*2 + 0];
                ykVec[i] = FT_Cr[i*2 + 1];
            }            
            break;
        
        case 2: // Inverse of Structure factor
            for (i=0; i<nrows; i++){
                rkVec[i] = StructFactor[i*2 + 0];
                ykVec[i] = 1.0 / StructFactor[i*2 + 1];
            }
            break;
        
        case 3: // Radial distribution function g(r)
            for (i=0; i<nrows; i++){
                rkVec[i] = Gr_data[i*2 + 0];
                ykVec[i] = Gr_data[i*2 + 1];
            }
            break;
        
        default : // Structure factor
            for (i=0; i<nrows; i++){
                rkVec[i] = StructFactor[i*2 + 0];
                ykVec[i] = StructFactor[i*2 + 1];
            }
            break;
    }

    free(r);
    free(q);
    free(U);
    free(Up);
    free(sigmaVec);
    free(StructFactor);
    free(FT_Cr);
    free(Gr_data);
    free(folderName);

    return 0;
}

/**
 * @brief Interpolates data using GSL Steffen splines.
 *
 * @param xInput Input x values.
 * @param yInput Input y values.
 * @param xOutput Output x values (interpolation points).
 * @param yOutput Output y values (interpolated results).
 * @param nrowsInput Number of input points.
 * @param nrowsOutput Number of output points.
 */
void interpolationFunc(double *xInput, double *yInput, double *xOutput, double *yOutput, int nrowsInput, int nrowsOutput) {
    
    size_t i;
    const size_t N = (size_t) nrowsInput;

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline_steffen = gsl_spline_alloc(gsl_interp_steffen, N);

    gsl_spline_init(spline_steffen, xInput, yInput, N);

    for (i = 0; i < nrowsOutput; ++i) {
        yOutput[i] = gsl_spline_eval(spline_steffen, xOutput[i], acc);
    }

    gsl_spline_free(spline_steffen);
    gsl_interp_accel_free(acc);

    return;
}

/**
 * @brief Generates a unique folder ID based on the current timestamp.
 *
 * @return A string containing the formatted timestamp.
 */
char *getFolderID() {
    char *fullTime, *day, *month, *year, *hour, *min, *sec;

    time_t currentTime;
    currentTime = time(NULL);

    fullTime = ctime(&currentTime);

    month = strtok(fullTime, " ");
    month = strtok(NULL, " ");
    day = strtok(NULL, " ");
    hour = strtok(NULL, ":");
    min = strtok(NULL, ":");
    sec = strtok(NULL, " ");
    year = strtok(NULL, "\n");

    int resultLength = strlen(day) + strlen(month) + strlen(year) +
                       strlen(hour) + strlen(min) + strlen(sec) + 6;

    char *fullTimeResult = (char *)malloc(resultLength);

    if (fullTimeResult == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    fullTimeResult[0] = '\0';
    strcat(fullTimeResult, day);
    strcat(fullTimeResult, month);
    strcat(fullTimeResult, year);
    strcat(fullTimeResult, "_");
    strcat(fullTimeResult, hour);
    strcat(fullTimeResult, min);
    strcat(fullTimeResult, sec);

    return fullTimeResult;
}

/**
 * @brief Checks if a directory exists.
 *
 * @param path Path to the directory.
 * @return true if the directory exists, false otherwise.
 */
bool directoryExists(char *path) {

    struct stat statbuf;
  
    if (stat(path, &statbuf) != -1) {
        return S_ISDIR(statbuf.st_mode);
    }
  
    return false;
}
