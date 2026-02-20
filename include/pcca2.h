#ifndef PCCA2_H
#define PCCA2_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gmp.h>

// Polynomial structure
// coeffs[i] represents the coefficient of x^i
typedef struct {
    mpq_t* coeffs;   // Coefficient array, size = degree+1, coeffs[i] is x^i coefficient
    int degree;

    // optional: store roots for verification
    mpq_t* roots;
    int nroots;
} Polynomial;
//-------------------------------------------------------
// methode util
// Create a polynomial
Polynomial* create_polynomial(int degree);

// Free polynomial memory
void free_polynomial(Polynomial *p);

// Free mpq_t*
void free_bound(mpq_t *bound);

// Copy(deep) a polynomial
Polynomial* copy_polynomial(Polynomial *p);

// Generate a random polynomial with coefficients in [min, max]
Polynomial* generate_random_polynomial(int n, long a, long b);

// afficher the result of the sturm
void afficher_bound(mpq_t* bound);

//-------------------------------------------------------
// polynomial operation
// PC=PA*PB
Polynomial* poly_mul_naive(Polynomial *A,Polynomial *B);

// Calculate the derivative of the polynomial
Polynomial* poly_derivative(Polynomial *p);

// recalculate the degree of poly
void poly_trim_degree(Polynomial *p);

// Calculate rem(A,B)
Polynomial* poly_remainder(Polynomial *A, Polynomial *B);

// (-1)*P
void poly_negative(Polynomial *p);

// // Calculate P(value)
double poly_calculate_sign(Polynomial *p, mpq_t x);

//-------------------------------------------------------
// basic calculation methods for sturm
// Calculate the cauchy bound of the roots of a polynomial
void cauchy_bound(mpq_t result, const Polynomial *p);

// Calculate the number of sign changes in a sequence
int nb_sign_change(double *l, int n);

// Check if the obtained results are correct (a corresponding root can be found for each interval).
int verify_interval(const mpq_t *bound, const Polynomial *p);

//-------------------------------------------------------
// methods of sturm sequences
// Find the interval corresponding to the root by splitting the input interval recursively
mpq_t* bound_recu(Polynomial **sturmSuite, int nbSuite, mpq_t inf, mpq_t sup);

// Find the interval of real roots of a polynomial by sturm
mpq_t* sturm_naif(Polynomial *p);

#endif