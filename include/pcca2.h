#ifndef PCCA2_H
#define PCCA2_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define EPSILON 1e-63

// Polynomial structure
// coeffs[i] represents the coefficient of x^i
typedef struct {
    double *coeffs;   // Coefficient array, size = degree+1, coeffs[i] is x^i coefficient
    int degree;

    double *roots;    // store roots for verification
    int nroots;
} Polynomial;

// Create a polynomial
Polynomial* create_polynomial(int degree);

// Free polynomial memory
void free_polynomial(Polynomial *p);

// Copy a polynomial
Polynomial* copy_polynomial(const Polynomial *p);

// Generate a random polynomial with coefficients in [min, max]
Polynomial* generate_random_polynomial(int n, double a, double b);

// afficher the result of the sturm
void afficher_bound(const double* bound);

// Calculate the derivative of the polynomial
Polynomial* poly_derivative(Polynomial *p);

// Calculate rem(A,B)
Polynomial* poly_remainder(Polynomial *A,Polynomial *B);

void poly_negative(Polynomial *p);

// Calculate P(value)
double poly_calculate(Polynomial *p, double value);

double cauchy_bound(Polynomial*p);

int nb_sign_change(double*l,int n);

int verify_interval(double * bound, Polynomial * p);

// function for sturm naif
double* bound_recu(Polynomial**sturmSuite, int nbSuite,double inf, double sup);

double* sturm_naif(Polynomial* p);

#endif