#include "structures_nonspherical.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

ProjectionMatrix* create_projection_matrix(int n_projections, int n_points) {
    ProjectionMatrix *pm = malloc(sizeof(ProjectionMatrix));
    if (!pm) return NULL;

    pm->n_projections = n_projections;
    pm->n_points = n_points;

    // Allocate array of pointers
    pm->data = malloc(n_projections * sizeof(double*));
    pm->labels = malloc(n_projections * sizeof(char*));

    if (!pm->data || !pm->labels) {
        if (pm->data) free(pm->data);
        if (pm->labels) free(pm->labels);
        free(pm);
        return NULL;
    }

    // Allocate rows
    for (int i = 0; i < n_projections; i++) {
        pm->data[i] = calloc(n_points, sizeof(double)); // Initialize to 0
        pm->labels[i] = NULL;
        if (!pm->data[i]) {
            // Rollback
            for (int k = 0; k < i; k++) free(pm->data[k]);
            free(pm->data);
            free(pm->labels);
            free(pm);
            return NULL;
        }
    }

    return pm;
}

void free_projection_matrix(ProjectionMatrix *pm) {
    if (!pm) return;

    if (pm->data) {
        for (int i = 0; i < pm->n_projections; i++) {
            if (pm->data[i]) free(pm->data[i]);
        }
        free(pm->data);
    }

    if (pm->labels) {
        for (int i = 0; i < pm->n_projections; i++) {
            if (pm->labels[i]) free(pm->labels[i]);
        }
        free(pm->labels);
    }

    free(pm);
}

void set_projection_label(ProjectionMatrix *pm, int index, const char *label) {
    if (!pm || index < 0 || index >= pm->n_projections) return;

    if (pm->labels[index]) free(pm->labels[index]);
    
    pm->labels[index] = strdup(label);
}
