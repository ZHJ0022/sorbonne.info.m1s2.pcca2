// Create a polynomial
#include "../include/pcca2.h"

Polynomial* create_polynomial(int degree) {
    Polynomial *p = (Polynomial*)malloc(sizeof(Polynomial));
    if (!p) return NULL;

    p->degree = degree;
    p->coeffs = (double*)calloc(degree + 1, sizeof(double));
    if (!p->coeffs) {
        free(p);
        return NULL;
    }

    return p;
}

// Free polynomial memory
void free_polynomial(Polynomial *p) {
    if (p) {
        if (p->coeffs) {
            free(p->coeffs);
        }
        free(p);
    }
}

// Copy a polynomial
Polynomial* copy_polynomial(const Polynomial *p) {
    if (!p) return NULL;

    Polynomial *copy = create_polynomial(p->degree);
    if (!copy) return NULL;

    memcpy(copy->coeffs, p->coeffs, (p->degree + 1) * sizeof(double));
    return copy;
}

// Generate a random polynomial
Polynomial* generate_random_polynomial(int degree, double min, double max) {
    Polynomial *p = create_polynomial(degree);
    if (!p) return NULL;

    for (int i = 0; i <= degree; i++) {
        // Generate a random number in the range [min, max]
        double random = (double)rand() / RAND_MAX;
        p->coeffs[i] = min + random * (max - min);
    }

    return p;
}