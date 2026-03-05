#ifndef STRUCTURES_NONSPHERICAL_H
#define STRUCTURES_NONSPHERICAL_H

#include <gsl/gsl_vector.h>

/**
 * @brief Structure to hold Rotational Invariant Projections.
 * 
 * For Dipolar Hard Spheres (DHS), we typically need 3 projections:
 * Index 0: h000 (Spherical/Density part)
 * Index 1: h110 (Angle dependent part 1)
 * Index 2: h112 (Angle dependent part 2 - Dipolar interaction)
 */
typedef struct {
    int n_projections;      // Number of projections (e.g., 3 for DHS)
    int n_points;           // Number of spatial points (r or k)
    double **data;          // 2D array [n_projections][n_points]
    char **labels;          // Labels for projections (e.g., "000", "110", "112")
} ProjectionMatrix;

/**
 * @brief Allocates memory for a ProjectionMatrix.
 * 
 * @param n_projections Number of projections.
 * @param n_points Number of spatial/k points.
 * @return Pointer to the allocated ProjectionMatrix.
 */
ProjectionMatrix* create_projection_matrix(int n_projections, int n_points);

/**
 * @brief Frees memory of a ProjectionMatrix.
 * 
 * @param pm Pointer to the ProjectionMatrix to free.
 */
void free_projection_matrix(ProjectionMatrix *pm);

/**
 * @brief Sets a label for a specific projection index.
 */
void set_projection_label(ProjectionMatrix *pm, int index, const char *label);

#endif /* STRUCTURES_NONSPHERICAL_H */
