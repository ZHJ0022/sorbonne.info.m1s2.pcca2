#ifndef PCCA2_H
#define PCCA2_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gmp.h>

#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq_mat.h>

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
// methode util(util.c)
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
// polynomial operation (poly.c)
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

// Calculate P(value) with Hornor
double poly_calculate_sign(Polynomial *p, mpq_t x);

// Convert Polynomial to fmpq_poly_t
void polynomialConv(fmpq_poly_t poly, Polynomial *p);

// Calculate P(value) with FLINT
int poly_sign_flint(fmpq_poly_t poly, mpq_t x);

// 0-> poly non-zero, 1->zero
int poly_is_zero(const Polynomial *p);

// get the lcoef
void poly_get_lc(mpq_t out, const Polynomial *p);

// 0-> non monic, 1-> monic(lcoef=1)
int poly_is_monic(const Polynomial *p);

// all coef *c
void poly_scalar_mul(Polynomial *p, mpq_t c);

// all coef /c
int poly_scalar_div(Polynomial *p, const mpq_t c);

//-------------------------------------------------------
// basic calculation methods for sturm/sturm havicht (calculate.c)
// Calculate the cauchy bound of the roots of a polynomial
void cauchy_bound(mpq_t result, Polynomial *p);

// Calculate the number of sign changes in a sequence
int nb_sign_change(double *l, int n);

// Check if the obtained results are correct (a corresponding root can be found for each interval).
int verify_interval(mpq_t *bound, Polynomial *p);

// (-1)^(k(k+1)/2)
int delta_sign(int k);

//-------------------------------------------------------
// methods of sturm sequences
// Find the interval corresponding to the root by splitting the input interval recursively
mpq_t* bound_recu(Polynomial **sturmSuite, int nbSuite, mpq_t inf, mpq_t sup, int nbI);

// Find the interval of real roots of a polynomial by sturm
mpq_t* sturm_naif(Polynomial *p, int nbI);

//-------------------------------------------------------
// (subesultant.c)

// create a new matrix of mpq_t
mpq_t** alloc_matrix(int rows, int cols);

// free the matrix of mpq_t
void free_matrix(mpq_t **M, int rows, int cols);

// create Sylvester matrix of P and S
mpq_t** create_sylv_j_matrix(Polynomial *P, Polynomial *S, int j, int *out_rows, int *out_cols);

// convert mpq_t(GMP) to fmpq_t(FLINT)
static void mpq_to_fmpq(fmpq_t y, const mpq_t x);

// convert fmpq_t(FLINT) to mpq_t(GMP)
static void fmpq_to_mpq(mpq_t y, const fmpq_t x);

// calculate the determinants of sylv matirx
void sylv_j_det(mpq_t out, mpq_t **M, int rows, int cols, int d);

// calculate the coefficients of the poly by the determinants of sylv matrix
Polynomial* subresultant_j(Polynomial *P, Polynomial *S, int j);

#endif