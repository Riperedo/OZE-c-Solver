#include "facdes2Y.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>

// =========================================================
// Función para Desplegar las Opciones de Potencial
// =========================================================
void display_potential_options() {
    printf("\n");
    printf("=========================================================================\n");
    printf("                  CATÁLOGO DE OPCIONES DE POTENCIALES\n");
    printf("=========================================================================\n");
    printf(" Caso | Nombre del Potencial      | Ecuación o Forma Característica\n");
    printf("------|---------------------------|---------------------------------------\n");
    printf("  1   | INVERSE POWER LAW (IPL)   | U = T* (sigma/r)^(lambda)\n");
    printf("  2   | TRUNCATED LENNARD-JONES   | Repulsivo (LJ Truncado en el minimo)\n");
    printf("  3   | TRUNCATED LENNARD-JONES 2 | LJ con minimo en r=sigma\n");
    printf("  4   | DOUBLE YUKAWA             | Atractivo + Repulsivo\n");
    printf("  5   | ATRACTIVE YUKAWA          | U = -T* exp[-lambda (r-1)]/r\n");
    printf("  6   | REPULSIVE YUKAWA          | U = T* exp[-lambda (r-1)]/r\n");
    printf("  7   | HARD SPHERE (HS)          | U = 0 (r > sigma)\n");
    printf("  8   | SHOULDER FUNCTION         | Potencial tipo 'Hombro' (Step Potencial)\n");
    printf("  9   | DOWN-HILL FUNCTION        | Potencial lineal decreciente\n");
    printf("  10  | GAUSSIAN CORE MODEL       | U = T* exp(- (r/sigma)^2 )\n");
    printf("  11  | RAMP (STEP FUNCTION)      | U lineal decreciente (tipo rampa)\n");
    printf("  12  | STEP FUNCTION (Soft Core) | U = E * (1 - r/sigma)^n (similar a Hertzian, pero se asume Step Function)\n");
    printf("  13  | HERTZIAN POTENTIAL (n=2.5)| U = E * (1 - r/sigma)^2.5 (r < sigma)\n");
    printf("-------------------------------------------------------------------------\n");
    printf("\n");
    printf("Ejemplo de uso: ./facdes_solver --closure HNC --potential 13 ...\n\n");
}

// Definición de una función auxiliar para imprimir el uso correcto del programa
void print_usage(const char *prog_name) {
    fprintf(stderr, "\nUso: %s [OPCION] [VALOR] ...\n\n", prog_name);
    fprintf(stderr, "Calcula el factor de estructura S(k) para un sistema coloidal.\n\n");
    fprintf(stderr, "Opciones requeridas:\n");
    fprintf(stderr, "  --closure   <HNC|RY>       Cierre termodinámico (HNC o RY).\n");
    fprintf(stderr, "  --potential <int>          ID del potencial de interacción (e.g., 1, 2, ...).\n");
    fprintf(stderr, "  --volfactor <double>       Factor de volumen (e.g., 0.1, 0.5).\n");
    fprintf(stderr, "  --temp      <double>       Temperatura T (e.g., 1.0).\n");
    fprintf(stderr, "  --nodes     <int>          Número de nodos internos para el cálculo (nodesFacdes2Y).\n");
    fprintf(stderr, "  --knodes    <int>          Número de nodos para el vector k de salida (k->size).\n");
    fprintf(stderr, "\nOpciones opcionales:\n");
    fprintf(stderr, "  --temp2     <double>       Segunda temperatura T2 (e.g., 1.0, por defecto 1.0).\n");
    fprintf(stderr, "  --lambda_a  <double>       Parámetro lambda_a (e.g., 0.1, por defecto 0.0).\n");
    fprintf(stderr, "  --lambda_r  <double>       Parámetro lambda_r (e.g., 0.1, por defecto 0.0).\n");
    fprintf(stderr, "\nEjemplo:\n");
    fprintf(stderr, "  %s--closure HNC --potential 7 --volfactor 0.2 --temp 1.0 --nodes 2048 --knodes 1024\n\n", prog_name);
}

// Función para imprimir ayuda específica según el potencial seleccionado
void print_potential_help(int potentialID) {
    printf("\n");
    printf("=========================================================================\n");
    printf("                  AYUDA ESPECÍFICA PARA EL POTENCIAL ID: %d\n", potentialID);
    printf("=========================================================================\n");

    switch (potentialID) {
        case 1: // INVERSE POWER LAW (IPL)
            printf("Potencial: INVERSE POWER LAW (IPL)\n");
            printf("Ecuación: U = T* (sigma/r)^(lambda)\n");
            printf("Parámetros REQUERIDOS:\n");
            printf("  --volfactor <double> : Factor de volumen\n");
            printf("  --temp      <double> : Temperatura T*\n");
            printf("  --lambda_a  <double> : Exponente lambda (repulsivo)\n");
            printf("  --nodes     <int>    : Nodos espaciales\n");
            printf("  --knodes    <int>    : Nodos en espacio k\n");
            printf("Ejemplo: \n");
            printf("./facdes_solver --closure HNC --potential 1 --volfactor 0.2 --temp 1.0 --lambda_a 12.0 --nodes 4096 --knodes 1024\n");
            break;

        case 2: // TRUNCATED LENNARD-JONES
        case 3: // TRUNCATED LENNARD-JONES 2
            printf("Potencial: TRUNCATED LENNARD-JONES (Tipo 1 o 2)\n");
            printf("Ecuación: Lennard-Jones truncado/desplazado\n");
            printf("Parámetros REQUERIDOS:\n");
            printf("  --volfactor <double> : Factor de volumen\n");
            printf("  --temp      <double> : Temperatura T*\n");
            printf("  --nodes     <int>    : Nodos espaciales\n");
            printf("  --knodes    <int>    : Nodos en espacio k\n");
            printf("Ejemplo: \n");
            printf("./facdes_solver --closure HNC --potential 2 --volfactor 0.25 --temp 1.0 --nodes 4096 --knodes 1024\n");
            break;

        case 4: // DOUBLE YUKAWA
            printf("Potencial: DOUBLE YUKAWA\n");
            printf("Ecuación: Atractivo + Repulsivo\n");
            printf("Parámetros REQUERIDOS:\n");
            printf("  --volfactor <double> : Factor de volumen\n");
            printf("  --temp      <double> : Temperatura T* (Atractiva)\n");
            printf("  --temp2     <double> : Temperatura T2* (Repulsiva)\n");
            printf("  --lambda_a  <double> : Lambda atractivo\n");
            printf("  --lambda_r  <double> : Lambda repulsivo\n");
            printf("  --nodes     <int>    : Nodos espaciales\n");
            printf("  --knodes    <int>    : Nodos en espacio k\n");
            printf("Ejemplo: \n");
            printf("./facdes_solver --closure HNC --potential 4 --volfactor 0.2 --temp 1.0 --temp2 0.5 --lambda_a 1.8 --lambda_r 5.0 --nodes 4096 --knodes 1024\n");
            break;

        case 5: // ATRACTIVE YUKAWA
            printf("Potencial: YUKAWA Atractivo\n");
            printf("Ecuación: U ~ exp(-lambda r)/r\n");
            printf("Parámetros REQUERIDOS:\n");
            printf("  --volfactor <double> : Factor de volumen\n");
            printf("  --temp      <double> : Temperatura T*\n");
            printf("  --lambda_a  <double> : Lambda\n");
            printf("  --nodes     <int>    : Nodos espaciales\n");
            printf("  --knodes    <int>    : Nodos en espacio k\n");
            printf("Ejemplo: \n");
            printf("./facdes_solver --closure HNC --potential 5 --volfactor 0.2 --temp 1.0 --lambda_a 1.8 --nodes 4096 --knodes 1024\n");
        case 6: // REPULSIVE YUKAWA
            printf("Potencial: YUKAWA Repulsivo\n");
            printf("Ecuación: U ~ exp(-lambda r)/r\n");
            printf("Parámetros REQUERIDOS:\n");
            printf("  --volfactor <double> : Factor de volumen\n");
            printf("  --temp      <double> : Temperatura T*\n");
            printf("  --lambda_a  <double> : Lambda\n");
            printf("  --nodes     <int>    : Nodos espaciales\n");
            printf("  --knodes    <int>    : Nodos en espacio k\n");
            printf("Ejemplo: \n");
            printf("./facdes_solver --closure HNC --potential 6 --volfactor 0.2 --temp 1.0 --lambda_a 1.8 --nodes 4096 --knodes 1024\n");
            break;

        case 7: // HARD SPHERE (HS)
            printf("Potencial: HARD SPHERE (HS)\n");
            printf("Ecuación: U = infinito si r < sigma, 0 si r > sigma\n");
            printf("Parámetros REQUERIDOS:\n");
            printf("  --volfactor <double> : Factor de volumen\n");
            printf("  --nodes     <int>    : Nodos espaciales\n");
            printf("  --knodes    <int>    : Nodos en espacio k\n");
            printf("Nota: La temperatura no afecta la estructura de HS puros, pero se requiere un valor dummy (e.g. 1.0).\n");
            printf("Ejemplo: \n");
            printf("./facdes_solver --closure HNC --potential 7 --volfactor 0.4 --temp 1.0 --nodes 4096 --knodes 1024\n");
            break;

        case 8: // SHOULDER FUNCTION
            printf("Potencial: SHOULDER FUNCTION (Step Potential)\n");
            printf("Ecuación: U = T* lambda   para sigma < r < T2\n");
            printf("Parámetros REQUERIDOS:\n");
            printf("  --volfactor <double> : Factor de volumen\n");
            printf("  --temp      <double> : Temperatura T*\n");
            printf("  --lambda_a  <double> : Altura del escalón (lambda)\n");
            printf("  --temp2     <double> : Ancho del escalón (T2)\n");
            printf("  --nodes     <int>    : Nodos espaciales\n");
            printf("  --knodes    <int>    : Nodos en espacio k\n");
            printf("Ejemplo: \n");
            printf("./facdes_solver --closure HNC --potential 8 --volfactor 0.3 --temp 1.0 --lambda_a 0.5 --temp2 1.5 --nodes 4096 --knodes 1024\n");
            break;

        case 9: // DOWN-HILL FUNCTION
            printf("Potencial: DOWN-HILL FUNCTION\n");
            printf("Ecuación: U lineal decreciente\n");
            printf("Parámetros REQUERIDOS:\n");
            printf("  --volfactor <double> : Factor de volumen\n");
            printf("  --temp      <double> : Temperatura T*\n");
            printf("  --lambda_a  <double> : Altura (lambda)\n");
            printf("  --temp2     <double> : Ancho (T2)\n");
            printf("  --nodes     <int>    : Nodos espaciales\n");
            printf("  --knodes    <int>    : Nodos en espacio k\n");
            printf("Ejemplo: \n");
            printf("./facdes_solver --closure HNC --potential 9 --volfactor 0.3 --temp 1.0 --lambda_a 0.5 --temp2 1.5 --nodes 4096 --knodes 1024\n");
            break;

        case 10: // GAUSSIAN CORE MODEL
            printf("Potencial: GAUSSIAN CORE MODEL\n");
            printf("Ecuación: U = T* exp(- (r/sigma)^2 )\n");
            printf("Parámetros REQUERIDOS:\n");
            printf("  --volfactor <double> : Factor de volumen\n");
            printf("  --temp      <double> : Temperatura T*\n");
            printf("  --nodes     <int>    : Nodos espaciales\n");
            printf("  --knodes    <int>    : Nodos en espacio k\n");
            printf("Ejemplo: \n");
            printf("./facdes_solver --closure HNC --potential 10 --volfactor 0.5 --temp 1.0 --nodes 4096 --knodes 1024\n");
            break;

        case 11: // RAMP (STEP FUNCTION)
            printf("Potencial: RAMP (STEP FUNCTION)\n");
            printf("Ecuación: U lineal decreciente (tipo rampa)\n");
            printf("Parámetros REQUERIDOS:\n");
            printf("  --volfactor <double> : Factor de volumen\n");
            printf("  --temp      <double> : Temperatura T*\n");
            printf("  --nodes     <int>    : Nodos espaciales\n");
            printf("  --knodes    <int>    : Nodos en espacio k\n");
            printf("Ejemplo: \n");
            printf("./facdes_solver --closure HNC --potential 11 --volfactor 0.3 --temp 1.0 --nodes 4096 --knodes 1024\n");
            break;

        case 12: // STEP FUNCTION (Soft Core)
            printf("Potencial: STEP FUNCTION (Soft Core)\n");
            printf("Ecuación: U = E * (1 - r/sigma)^n\n");
            printf("Parámetros REQUERIDOS:\n");
            printf("  --volfactor <double> : Factor de volumen\n");
            printf("  --temp      <double> : Temperatura T* (E)\n");
            printf("  --lambda_a  <double> : Exponente n (lambda)\n");
            printf("  --temp2     <double> : Ancho (T2)\n");
            printf("  --nodes     <int>    : Nodos espaciales\n");
            printf("  --knodes    <int>    : Nodos en espacio k\n");
            printf("Ejemplo: \n");
            printf("./facdes_solver --closure HNC --potential 12 --volfactor 0.3 --temp 1.0 --lambda_a 2.0 --temp2 1.5 --nodes 4096 --knodes 1024\n");
            break;

        case 13: // HERTZIAN POTENTIAL
            printf("Potencial: HERTZIAN POTENTIAL (n=2.5)\n");
            printf("Ecuación: U = E * (1 - r/sigma)^2.5 (r < sigma)\n");
            printf("Parámetros REQUERIDOS:\n");
            printf("  --volfactor <double> : Factor de volumen\n");
            printf("  --temp      <double> : Pre-factor de energía E (o T*)\n");
            printf("  --nodes     <int>    : Nodos espaciales\n");
            printf("  --knodes    <int>    : Nodos en espacio k\n");
            printf("Ejemplo: \n");
            printf("./facdes_solver --closure HNC --potential 13 --volfactor 0.3 --temp 1.0 --nodes 4096 --knodes 1024\n");
            break;

        default:
            printf("Potencial ID %d no tiene ayuda específica detallada aún.\n", potentialID);
            printf("Revise la lista general de potenciales.\n");
            printf("Parámetros típicamente requeridos: --volfactor, --temp, --nodes, --knodes.\n");
            break;
    }
    printf("=========================================================================\n\n");
}

int main(int argc, char *argv[]) {
    // Si el usuario no proporciona argumentos, muestra el uso y las opciones.
    if (argc < 2) {
        printf("Uso: %s --closure [HNC/PY/RY] --potential [ID] --volfactor [...] ...\n", argv[0]);
        display_potential_options(); 
        return 1;
    }

    // Valores predeterminados y variables de entrada
    char *closure_str = NULL;
    int potentialNumber = -1;
    double volumeFactor = -1.0;
    double Temperature = -1.0;
    int nodesFacdes2Y = -1;
    int k_nodes = -1;
    
    // Valores opcionales
    double Temperature2 = 1.0;
    double lambda_a = 0.0;
    double lambda_r = 0.0;
    
    // Parseo de argumentos de línea de comandos
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
        } else if (strcmp(argv[i], "--nodes") == 0 && i + 1 < argc) {
            nodesFacdes2Y = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--knodes") == 0 && i + 1 < argc) {
            k_nodes = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return EXIT_SUCCESS;
        } else {
            fprintf(stderr, "Error: Opción desconocida o incompleta: %s\n", argv[i]);
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    // Validación de argumentos requeridos
    if (closure_str == NULL || potentialNumber == -1 || volumeFactor == -1.0 || 
        Temperature == -1.0 || nodesFacdes2Y == -1 || k_nodes == -1) {
        
        // Si tenemos un ID de potencial pero faltan otros argumentos, mostramos la ayuda específica
        if (potentialNumber != -1) {
            print_potential_help(potentialNumber);
            return EXIT_FAILURE;
        }

        fprintf(stderr, "Error: Faltan argumentos requeridos.\n");
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    
    if (strcmp(closure_str, "HNC") != 0 && strcmp(closure_str, "RY") != 0) {
        fprintf(stderr, "Error: Cierre ('%s') no válido. Use 'HNC' o 'RY'.\n", closure_str);
        return EXIT_FAILURE;
    }

    // Preparación del vector k (espacio de Fourier)
    gsl_vector *k_vec = gsl_vector_alloc(k_nodes);
    if (k_vec == NULL) {
        fprintf(stderr, "Error de asignación de memoria para gsl_vector k.\n");
        return EXIT_FAILURE;
    }
    
    // Preparación del vector r (espacio real)
    gsl_vector *r_vec = gsl_vector_alloc(k_nodes);
    if (r_vec == NULL) {
        fprintf(stderr, "Error de asignación de memoria para gsl_vector r.\n");
        gsl_vector_free(k_vec);
        return EXIT_FAILURE;
    }
    
    // Inicialización simple de k: valores espaciados uniformemente
    // Desde un k_min pequeño (para evitar la singularidad en 0) hasta un k_max.
    double k_max = 10.0; 
    double k_min = k_max / (double)k_nodes; // Valor pequeño > 0
    double dk = (k_max - k_min) / (double)(k_nodes - 1);
    
    for (int i = 0; i < k_nodes; i++) {
        gsl_vector_set(k_vec, i, k_min + i * dk);
    }
    
    // Inicialización del vector r (espacio real)
    double r_min = 0.01;
    double r_max = 10.0;
    double dr = (r_max - r_min) / (double)(k_nodes - 1);
    
    for (int i = 0; i < k_nodes; i++) {
        gsl_vector_set(r_vec, i, r_min + i * dr);
    }
    
    // Vector de salida para S(k)
    double *sk_output = malloc(k_nodes * sizeof(double));
    if (sk_output == NULL) {
        fprintf(stderr, "Error de asignación de memoria para sk_output.\n");
        gsl_vector_free(k_vec);
        gsl_vector_free(r_vec);
        return EXIT_FAILURE;
    }
    
    // Vector de salida para g(r)
    double *gr_output = malloc(k_nodes * sizeof(double));
    if (gr_output == NULL) {
        fprintf(stderr, "Error de asignación de memoria para gr_output.\n");
        gsl_vector_free(k_vec);
        gsl_vector_free(r_vec);
        free(sk_output);
        return EXIT_FAILURE;
    }

    printf("Iniciando cálculo...\n");
    printf("Cierre: %s, Potencial: %d, phi: %.4f, T: %.4f, N_calc: %d, N_k: %d\n", 
           closure_str, potentialNumber, volumeFactor, Temperature, nodesFacdes2Y, k_nodes);

    // Llamada a las funciones correspondientes
    if (strcmp(closure_str, "HNC") == 0) {
        printf("\n# Calculando S(k)...\n");
        sk_HNC(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, 
               k_vec, sk_output, potentialNumber, nodesFacdes2Y);
        
        printf("\n# Calculando g(r)...\n");
        gr_HNC(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, 
               r_vec, gr_output, potentialNumber, nodesFacdes2Y);
    } else { // Cierre "RY"
        printf("\n# Calculando S(k)...\n");
        sk_RY(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, 
              k_vec, sk_output, potentialNumber, nodesFacdes2Y);
        
        printf("\n# Calculando g(r)...\n");
        gr_RY(volumeFactor, Temperature, Temperature2, lambda_a, lambda_r, 
              r_vec, gr_output, potentialNumber, nodesFacdes2Y);
    }

    
    // Liberar memoria
    free(sk_output);
    free(gr_output);
    gsl_vector_free(k_vec);
    gsl_vector_free(r_vec);
    
    return EXIT_SUCCESS;
}

// compilar
// gcc -Wall -o facdes_solver main.c facdes2Y.c math_aux.c structures.c -lgsl -lgslcblas -lm