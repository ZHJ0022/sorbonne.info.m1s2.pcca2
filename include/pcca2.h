#ifndef PCCA2_H
#define PCCA2_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Polynomial structure
// coeffs[i] represents the coefficient of x^i
typedef struct {
    double *coeffs;    // Coefficient array
    int degree;        // Degree of the polynomial
} Polynomial;

// Create a polynomial
Polynomial* create_polynomial(int degree);

// Free polynomial memory
void free_polynomial(Polynomial *p);

// Copy a polynomial
Polynomial* copy_polynomial(const Polynomial *p);

// Generate a random polynomial with coefficients in [min, max]
Polynomial* generate_random_polynomial(int degree, double min, double max);

#endif