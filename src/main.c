#include "facdes2Y.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>

// Forward declaration of the new solver
void solver_dipolar(int closureID, double temp, double rho, double dipole_moment, 
                    int nodes, double rmax, const char *output_dir,
                    double temp_start, int ramp_steps, const char *init_sk_file);
void solver_mode2_core(int closureID, double temp, double rho, double dipole_moment, 
                   int nodes, double rmax, const char *output_dir);

// =========================================================
// Function to Display the Catalog of Potential Options
// =========================================================
void display_potential_options() {
    printf("\n");
    printf("=========================================================================\n");
    printf("                  POTENTIAL INTERACTION CATALOG\n");
    printf("=========================================================================\n");
    printf(" Case | Potential Name            | Characteristic Equation / Form\n");
    printf("------|---------------------------|---------------------------------------\n");
    printf("  1   | INVERSE POWER LAW (IPL)   | U = T* (sigma/r)^(lambda)\n");
    printf("  2   | TRUNCATED LENNARD-JONES   | Repulsive (LJ truncated at minimum)\n");
    printf("  3   | TRUNCATED LENNARD-JONES 2 | LJ with minimum shifted to r = sigma\n");
    printf("  4   | DOUBLE YUKAWA             | Attractive + Repulsive\n");
    printf("  5   | ATTRACTIVE YUKAWA         | U = -T* exp[-lambda (r-1)]/r\n");
    printf("  6   | REPULSIVE YUKAWA          | U = T* exp[-lambda (r-1)]/r\n");
    printf("  7   | HARD SPHERE (HS)          | U = infinity (r < sigma), 0 (r > sigma)\n");
    printf("  8   | SHOULDER FUNCTION         | Square shoulder / step potential\n");
    printf("  9   | DOWN-HILL FUNCTION        | Linearly decreasing repulsive potential\n");
    printf("  10  | GAUSSIAN CORE MODEL       | U = T* exp(- (r/sigma)^2 )\n");
    printf("  11  | RAMP (STEP FUNCTION)      | Linearly decreasing ramp potential\n");
    printf("  12  | STEP FUNCTION (Soft Core) | U = E * (1 - r/sigma)^n\n");
    printf("  13  | HERTZIAN POTENTIAL (n=2.5)| U = E * (1 - r/sigma)^2.5 (r < sigma)\n");
    printf("  14  | DIPOLAR HARD SPHERES      | Hard Spheres + Point Dipole (MSA/LHNC/QHNC/RHNC)\n");
    printf("  15  | DHS EXTENDED (Modes 2)    | Coupled higher rotational invariants (m,n <= 2)\n");
    printf("  16  | SOFT SHOULDER POTENTIAL   | U = 0.5 * E * (1 - tanh(alpha * (r - lambda)))\n");
    printf("-------------------------------------------------------------------------\n");
    printf("\n");
    printf("Example usage: ./build/facdes_solver --closure HNC --potential 13 ...\n\n");
}

// Definition of helper function to print usage instructions
void print_usage(const char *prog_name) {
    fprintf(stderr, "\nUsage: %s [OPTION] [VALUE] ...\n\n", prog_name);
    fprintf(stderr, "Solves the Ornstein-Zernike Equation to calculate structure factor S(k) and radial distribution g(r).\n\n");
    fprintf(stderr, "Required options:\n");
    fprintf(stderr, "  --closure   <HNC|RY|PY|MSA|LHNC|QHNC|RHNC>  Thermodynamic closure relation.\n");
    fprintf(stderr, "  --potential <int>                           Interaction potential ID (1 to 16).\n");
    fprintf(stderr, "  --volfactor <double>                        Packing fraction / volume factor eta (e.g. 0.1, 0.4).\n");
    fprintf(stderr, "  --temp      <double>                        Reduced temperature T* (e.g. 1.0).\n");
    fprintf(stderr, "  --nodes     <int>                           Internal radial discretization nodes (e.g. 2048, 4096).\n");
    fprintf(stderr, "  --knodes    <int>                           Number of Fourier k-space nodes for isotropic output.\n");
    fprintf(stderr, "\nOptional parameters:\n");
    fprintf(stderr, "  --temp2     <double>       Secondary temperature T2* (default: 1.0).\n");
    fprintf(stderr, "  --lambda_a  <double>       Attractive range parameter (default: 0.0).\n");
    fprintf(stderr, "  --lambda_r  <double>       Repulsive range parameter (default: 0.0).\n");
    fprintf(stderr, "  --dipole    <double>       Reduced dipole moment mu* (required for potentials 14 and 15).\n");
    fprintf(stderr, "  --rmax      <double>       Maximum radial domain cutoff r_max (default: 15.0 for dipoles / 10.0 for spherical).\n");
    fprintf(stderr, "  --temp-start <double>      Initial high temperature for continuation annealing ramp.\n");
    fprintf(stderr, "  --temp-steps <int>         Number of intermediate continuation stages (default: 10).\n");
    fprintf(stderr, "  --ramp                     Enable automated geometric temperature continuation ramp.\n");
    fprintf(stderr, "  --init-sk   <string>       Input .dat file with analytical/prior structure factors (k, S000, S110, S112) for warm-start.\n");
    fprintf(stderr, "\nExamples:\n");
    fprintf(stderr, "  Spherical (Hertzian):  %s --closure HNC --potential 13 --volfactor 0.3 --temp 1.0 --nodes 4096 --knodes 1024\n", prog_name);
    fprintf(stderr, "  Dipolar (RHNC):        %s --closure RHNC --potential 14 --volfactor 0.418879 --temp 1.0 --dipole 1.6583 --nodes 4096 --rmax 30.0\n\n", prog_name);
}

// Helper function to print detailed parameter guidance per potential ID
void print_potential_help(int potentialID) {
    printf("\n");
    printf("=========================================================================\n");
    printf("                  DETAILED HELP FOR POTENTIAL ID: %d\n", potentialID);
    printf("=========================================================================\n");

    switch (potentialID) {
        case 1: // INVERSE POWER LAW (IPL)
            printf("Potential: INVERSE POWER LAW (IPL)\n");
            printf("Equation: U = T* (sigma/r)^(lambda)\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --temp      <double> : Temperature T*\n");
            printf("  --lambda_a  <double> : Repulsive exponent lambda\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("  --knodes    <int>    : k-space nodes\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure HNC --potential 1 --volfactor 0.2 --temp 1.0 --lambda_a 12.0 --nodes 4096 --knodes 1024\n");
            break;

        case 2: // TRUNCATED LENNARD-JONES
        case 3: // TRUNCATED LENNARD-JONES 2
            printf("Potential: TRUNCATED LENNARD-JONES (Type 1 or 2)\n");
            printf("Equation: Truncated/shifted Lennard-Jones\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --temp      <double> : Temperature T*\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("  --knodes    <int>    : k-space nodes\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure HNC --potential 2 --volfactor 0.25 --temp 1.0 --nodes 4096 --knodes 1024\n");
            break;

        case 4: // DOUBLE YUKAWA
            printf("Potential: DOUBLE YUKAWA\n");
            printf("Equation: Attractive + Repulsive Yukawa\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --temp      <double> : Attractive temperature T*\n");
            printf("  --temp2     <double> : Repulsive temperature T2*\n");
            printf("  --lambda_a  <double> : Attractive screening parameter\n");
            printf("  --lambda_r  <double> : Repulsive screening parameter\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("  --knodes    <int>    : k-space nodes\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure HNC --potential 4 --volfactor 0.2 --temp 1.0 --temp2 0.5 --lambda_a 1.8 --lambda_r 5.0 --nodes 4096 --knodes 1024\n");
            break;

        case 5: // ATTRACTIVE YUKAWA
            printf("Potential: ATTRACTIVE YUKAWA\n");
            printf("Equation: U ~ exp(-lambda r)/r\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --temp      <double> : Temperature T*\n");
            printf("  --lambda_a  <double> : Screening parameter lambda\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("  --knodes    <int>    : k-space nodes\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure HNC --potential 5 --volfactor 0.2 --temp 1.0 --lambda_a 1.8 --nodes 4096 --knodes 1024\n");
            break;
        case 6: // REPULSIVE YUKAWA
            printf("Potential: REPULSIVE YUKAWA\n");
            printf("Equation: U ~ exp(-lambda r)/r\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --temp      <double> : Temperature T*\n");
            printf("  --lambda_a  <double> : Screening parameter lambda\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("  --knodes    <int>    : k-space nodes\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure HNC --potential 6 --volfactor 0.2 --temp 1.0 --lambda_a 1.8 --nodes 4096 --knodes 1024\n");
            break;

        case 7: // HARD SPHERE (HS)
            printf("Potential: HARD SPHERE (HS)\n");
            printf("Equation: U = infinity (r < sigma), 0 (r > sigma)\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("  --knodes    <int>    : k-space nodes\n");
            printf("Note: Temperature does not alter pure HS structure, but a dummy value (e.g. 1.0) is passed.\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure HNC --potential 7 --volfactor 0.4 --temp 1.0 --nodes 4096 --knodes 1024\n");
            break;

        case 8: // SHOULDER FUNCTION
            printf("Potential: SHOULDER FUNCTION (Step Potential)\n");
            printf("Equation: U = T* lambda for sigma < r < T2\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --temp      <double> : Temperature T*\n");
            printf("  --lambda_a  <double> : Step height (lambda)\n");
            printf("  --temp2     <double> : Step width (T2)\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("  --knodes    <int>    : k-space nodes\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure HNC --potential 8 --volfactor 0.3 --temp 1.0 --lambda_a 0.5 --temp2 1.5 --nodes 4096 --knodes 1024\n");
            break;

        case 9: // DOWN-HILL FUNCTION
            printf("Potential: DOWN-HILL FUNCTION\n");
            printf("Equation: Linearly decreasing repulsive potential\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --temp      <double> : Temperature T*\n");
            printf("  --lambda_a  <double> : Height (lambda)\n");
            printf("  --temp2     <double> : Width (T2)\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("  --knodes    <int>    : k-space nodes\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure HNC --potential 9 --volfactor 0.3 --temp 1.0 --lambda_a 0.5 --temp2 1.5 --nodes 4096 --knodes 1024\n");
            break;

        case 10: // GAUSSIAN CORE MODEL
            printf("Potential: GAUSSIAN CORE MODEL\n");
            printf("Equation: U = T* exp(- (r/sigma)^2 )\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --temp      <double> : Temperature T*\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("  --knodes    <int>    : k-space nodes\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure HNC --potential 10 --volfactor 0.5 --temp 1.0 --nodes 4096 --knodes 1024\n");
            break;

        case 11: // RAMP (STEP FUNCTION)
            printf("Potential: RAMP (STEP FUNCTION)\n");
            printf("Equation: Linearly decreasing ramp\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --temp      <double> : Temperature T*\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("  --knodes    <int>    : k-space nodes\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure HNC --potential 11 --volfactor 0.3 --temp 1.0 --nodes 4096 --knodes 1024\n");
            break;

        case 12: // STEP FUNCTION (Soft Core)
            printf("Potential: STEP FUNCTION (Soft Core)\n");
            printf("Equation: U = E * (1 - r/sigma)^n\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --temp      <double> : Temperature T* (E)\n");
            printf("  --lambda_a  <double> : Exponent n\n");
            printf("  --temp2     <double> : Width (T2)\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("  --knodes    <int>    : k-space nodes\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure HNC --potential 12 --volfactor 0.3 --temp 1.0 --lambda_a 2.0 --temp2 1.5 --nodes 4096 --knodes 1024\n");
            break;

        case 13: // HERTZIAN POTENTIAL
            printf("Potential: HERTZIAN POTENTIAL (n=2.5)\n");
            printf("Equation: U = E * (1 - r/sigma)^2.5 (r < sigma)\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --temp      <double> : Energy prefactor E (or T*)\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("  --knodes    <int>    : k-space nodes\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure HNC --potential 13 --volfactor 0.3 --temp 1.0 --nodes 4096 --knodes 1024\n");
            break;

        case 14: // DIPOLAR HARD SPHERES
            printf("Potential: DIPOLAR HARD SPHERES\n");
            printf("Equation: Hard Spheres + Point Dipole Interaction\n");
            printf("Supported closures: MSA, LHNC, QHNC, RHNC\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction eta = (pi/6)*rho*sigma^3\n");
            printf("  --temp      <double> : Reduced temperature T*\n");
            printf("  --dipole    <double> : Dipole moment mu*\n");
            printf("  --nodes     <int>    : Real-space nodes (e.g. 4096)\n");
            printf("Optional flags: --rmax, --ramp, --temp-start, --temp-steps, --init-sk\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure RHNC --potential 14 --volfactor 0.418879 --temp 1.0 --dipole 1.6583 --nodes 4096 --rmax 30.0\n");
            break;

        case 15: // DHS EXTENDED (Modes 2)
            printf("Potential: DHS EXTENDED (Modes 2)\n");
            printf("Equation: Coupled higher rotational invariants (m,n <= 2)\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --temp      <double> : Temperature T*\n");
            printf("  --dipole    <double> : Dipole moment mu*\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure LHNC --potential 15 --volfactor 0.3 --temp 1.0 --dipole 1.5 --nodes 2048\n");
            break;

        case 16: // SOFT SHOULDER POTENTIAL
            printf("Potential: SOFT SHOULDER POTENTIAL\n");
            printf("Equation: U = 0.5 * E * (1 - tanh(alpha * (r - lambda)))\n");
            printf("Required arguments:\n");
            printf("  --volfactor <double> : Packing fraction\n");
            printf("  --temp      <double> : Energy scale E (or T*)\n");
            printf("  --lambda_a  <double> : Potential range (lambda)\n");
            printf("  --lambda_r  <double> : Shoulder steepness (alpha)\n");
            printf("  --nodes     <int>    : Real-space nodes\n");
            printf("  --knodes    <int>    : k-space nodes\n");
            printf("Example:\n");
            printf("./build/facdes_solver --closure HNC --potential 16 --volfactor 0.2 --temp 1.0 --lambda_a 2.0 --lambda_r 5.0 --nodes 2048 --knodes 1024\n");
            break;

        default:
            printf("Potential ID %d does not have detailed help documentation yet.\n", potentialID);
            printf("Please review the general potential catalog.\n");
            printf("Typically required parameters: --volfactor, --temp, --nodes, --knodes.\n");
            break;
    }
    printf("=========================================================================\n\n");
}

int main(int argc, char *argv[]) {
    // If user provides no arguments, show usage and options.
    if (argc < 2) {
        printf("Usage: %s --closure [HNC|RY|PY|MSA|LHNC|QHNC|RHNC] --potential [ID] --volfactor [...] ...\n", argv[0]);
        display_potential_options(); 
        return 1;
    }

    // Default values and input variables
    char *closure_str = NULL;
    int potentialNumber = -1;
    double volumeFactor = -1.0;
    double Temperature = -1.0;
    int nodesFacdes2Y = -1;
    int k_nodes = -1;
    
    // Optional values
    double Temperature2 = 1.0;
    double lambda_a = 0.0;
    double lambda_r = 0.0;
    double dipole_moment = 0.0;
    double rmax_val = 15.0; // Default rmax (15.0 for dipoles)
    double temp_start = -1.0; // Default: no temperature ramping
    int temp_steps = 10;      // Number of continuation stages
    int use_ramp = 0;         // Flag for automatic temperature ramping
    const char *init_sk_file = NULL; // Analytical structure factor input file
    
    // Parse command-line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--closure") == 0 && i + 1 < argc) {
            closure_str = argv[++i];
        } else if (strcmp(argv[i], "--potential") == 0 && i + 1 < argc) {
            potentialNumber = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--volfactor") == 0 && i + 1 < argc) {
            volumeFactor = atof(argv[++i]);
        } else if (strcmp(argv[i], "--temp") == 0 && i + 1 < argc) {
            Temperature = atof(argv[++i]);
        } else if (strcmp(argv[i], "--temp2") == 0 && i + 1 < argc) {
            Temperature2 = atof(argv[++i]);
        } else if (strcmp(argv[i], "--lambda_a") == 0 && i + 1 < argc) {
            lambda_a = atof(argv[++i]);
        } else if (strcmp(argv[i], "--lambda_r") == 0 && i + 1 < argc) {
            lambda_r = atof(argv[++i]);
        } else if (strcmp(argv[i], "--dipole") == 0 && i + 1 < argc) {
            dipole_moment = atof(argv[++i]);
        } else if (strcmp(argv[i], "--rmax") == 0 && i + 1 < argc) {
            rmax_val = atof(argv[++i]);
        } else if (strcmp(argv[i], "--temp-start") == 0 && i + 1 < argc) {
            temp_start = atof(argv[++i]);
        } else if (strcmp(argv[i], "--temp-steps") == 0 && i + 1 < argc) {
            temp_steps = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--ramp") == 0) {
            use_ramp = 1;
        } else if (strcmp(argv[i], "--init-sk") == 0 && i + 1 < argc) {
            init_sk_file = argv[++i];
        } else if (strcmp(argv[i], "--nodes") == 0 && i + 1 < argc) {
            nodesFacdes2Y = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--knodes") == 0 && i + 1 < argc) {
            k_nodes = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return EXIT_SUCCESS;
        } else {
            fprintf(stderr, "Error: Unknown or incomplete option: %s\n", argv[i]);
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    // For dipolar systems (potential 14 and 15), knodes defaults to nodesFacdes2Y
    if ((potentialNumber == 14 || potentialNumber == 15) && k_nodes == -1) {
        k_nodes = nodesFacdes2Y;
    }

    // Validation of required arguments
    if (closure_str == NULL || potentialNumber == -1 || volumeFactor == -1.0 || 
        Temperature == -1.0 || nodesFacdes2Y == -1 || k_nodes == -1) {
        
        // If potential ID is provided but other arguments are missing, show specific help
        if (potentialNumber != -1) {
            print_potential_help(potentialNumber);
            return EXIT_FAILURE;
        }

        fprintf(stderr, "Error: Missing required arguments.\n");
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    
    // Check for Dipolar Solver (Potential 14)
    if (potentialNumber == 14) {
        if (dipole_moment <= 0.0) {
            fprintf(stderr, "Error: For potential 14, --dipole must be > 0.\n");
            return EXIT_FAILURE;
        }

        int closure_id_int = 0; // Default MSA
        if (strcmp(closure_str, "LHNC") == 0 || strcmp(closure_str, "1") == 0) closure_id_int = 1;
        else if (strcmp(closure_str, "QHNC") == 0 || strcmp(closure_str, "2") == 0) closure_id_int = 2;
        else if (strcmp(closure_str, "RHNC") == 0 || strcmp(closure_str, "3") == 0) closure_id_int = 3;
        else if (strcmp(closure_str, "HNC") == 0) closure_id_int = 1; // Fallback to LHNC for "HNC" input
        else if (strcmp(closure_str, "MSA") == 0 || strcmp(closure_str, "0") == 0) closure_id_int = 0;
        else if (strcmp(closure_str, "RY") == 0) {
             fprintf(stderr, "Warning: RY is not implemented for Dipolar fluids yet. Using LHNC.\n");
             closure_id_int = 1;
        }

        // Handle temperature ramping defaults
        if (use_ramp && temp_start <= 0.0) {
            temp_start = (Temperature < 2.0) ? 5.0 : Temperature * 2.0;
        }

        // The input volumeFactor is the packing fraction eta. 
        // For hard spheres, eta = (pi/6) * rho * sigma^3. With sigma=1, rho = 6 * eta / pi.
        double rho = 6.0 * volumeFactor / M_PI;

        // Call the non-spherical dipolar solver with temperature continuation and analytical S(k) init support
        solver_dipolar(closure_id_int, Temperature, rho, dipole_moment, nodesFacdes2Y, rmax_val, "output", temp_start, temp_steps, init_sk_file);
        return EXIT_SUCCESS;
    }

    // Check for Extended Modes 2 Solver (Potential 15)
    if (potentialNumber == 15) {
        if (dipole_moment <= 0.0) {
            fprintf(stderr, "Error: For potential 15, --dipole must be > 0.\n");
            return EXIT_FAILURE;
        }

        int closure_id_int = 0; // Default MSA
        if (strcmp(closure_str, "LHNC") == 0) closure_id_int = 1;
        else if (strcmp(closure_str, "QHNC") == 0) closure_id_int = 2;
        else if (strcmp(closure_str, "HNC") == 0) closure_id_int = 1; 
        else if (strcmp(closure_str, "MSA") == 0) closure_id_int = 0;
        
        double rho = 6.0 * volumeFactor / M_PI;
        solver_mode2_core(closure_id_int, Temperature, rho, dipole_moment, nodesFacdes2Y, rmax_val, "output");
        return EXIT_SUCCESS;
    }

    if (strcmp(closure_str, "HNC") != 0 && strcmp(closure_str, "RY") != 0 && strcmp(closure_str, "PY") != 0) {
        fprintf(stderr, "Error: Invalid closure ('%s'). Use 'HNC', 'RY', or 'PY'.\n", closure_str);
        return EXIT_FAILURE;
    }

    // Allocation of k-space vector (Fourier space)
    gsl_vector *k_vec = gsl_vector_alloc(k_nodes);
    if (k_vec == NULL) {
        fprintf(stderr, "Memory allocation error for gsl_vector k.\n");
        return EXIT_FAILURE;
    }
    
    // Allocation of real-space r vector
    gsl_vector *r_vec = gsl_vector_alloc(k_nodes);
    if (r_vec == NULL) {
        fprintf(stderr, "Memory allocation error for gsl_vector r.\n");
        gsl_vector_free(k_vec);
        return EXIT_FAILURE;
    }
    
    // Linear discretization of k grid
    double k_max = 10.0; 
    double k_min = k_max / (double)k_nodes;
    double dk = (k_max - k_min) / (double)(k_nodes - 1);
    
    for (int i = 0; i < k_nodes; i++) {
        gsl_vector_set(k_vec, i, k_min + i * dk);
    }
    
    // Linear discretization of r grid
    double r_min = 0.01;
    double r_max = 10.0;
    double dr = (r_max - r_min) / (double)(k_nodes - 1);
    
    for (int i = 0; i < k_nodes; i++) {
        gsl_vector_set(r_vec, i, r_min + i * dr);
    }
    
    // Output array for S(k)
    double *sk_output = malloc(k_nodes * sizeof(double));
    if (sk_output == NULL) {
        fprintf(stderr, "Memory allocation error for sk_output.\n");
        gsl_vector_free(k_vec);
        gsl_vector_free(r_vec);
        return EXIT_FAILURE;
    }
    
    // Output array for g(r)
    double *gr_output = malloc(k_nodes * sizeof(double));
    if (gr_output == NULL) {
        fprintf(stderr, "Memory allocation error for gr_output.\n");
        gsl_vector_free(k_vec);
        gsl_vector_free(r_vec);
        free(sk_output);
        return EXIT_FAILURE;
    }

    printf("Starting calculation...\n");
    printf("Closure: %s, Potential: %d, phi: %.4f, T*: %.4f, N_nodes: %d, N_k: %d\n", 
           closure_str, potentialNumber, volumeFactor, Temperature, nodesFacdes2Y, k_nodes);

    // Dispatch to corresponding calculation functions
    if (strcmp(closure_str, "HNC") == 0) {
        printf("\n# Calculating S(k)...\n");
        sk_HNC(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, 
               k_vec, sk_output, potentialNumber, nodesFacdes2Y);
        
        printf("\n# Calculating g(r)...\n");
        gr_HNC(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, 
               r_vec, gr_output, potentialNumber, nodesFacdes2Y);
    } else if (strcmp(closure_str, "PY") == 0) {
        printf("\n# Calculating S(k)...\n");
        sk_PY(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, 
               k_vec, sk_output, potentialNumber, nodesFacdes2Y);
        
        printf("\n# Calculating g(r)...\n");
        gr_PY(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, 
               r_vec, gr_output, potentialNumber, nodesFacdes2Y);
    } else { // "RY" Closure
        printf("\n# Calculating S(k)...\n");
        sk_RY(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, 
               k_vec, sk_output, potentialNumber, nodesFacdes2Y);
        
        printf("\n# Calculating g(r)...\n");
        gr_RY(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, 
               r_vec, gr_output, potentialNumber, nodesFacdes2Y);
    }

    // Free memory
    free(sk_output);
    free(gr_output);
    gsl_vector_free(k_vec);
    gsl_vector_free(r_vec);
    
    return EXIT_SUCCESS;
}

// compilar
// gcc -Wall -o facdes_solver main.c facdes2Y.c math_aux.c structures.c -lgsl -lgslcblas -lm