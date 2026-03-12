#include "../include/pcca2.h"

// PC=PA*PB
Polynomial* poly_mul_naive(Polynomial *A,Polynomial *B) {
    if (!A || !B) return NULL;

    Polynomial *C = create_polynomial(A->degree + B->degree);
    if (!C) return NULL;

    mpq_t tmp;
    mpq_init(tmp);

    for (int i = 0; i <= A->degree; i++) {
        for (int j = 0; j <= B->degree; j++) {
            // tmp = A[i] * B[j]
            mpq_mul(tmp, A->coeffs[i], B->coeffs[j]);

            // C[i+j] += tmp
            mpq_add(C->coeffs[i + j], C->coeffs[i + j], tmp);
        }
    }

    mpq_clear(tmp);
    return C;
}

// Calculate the derivative of the polynomial
Polynomial* poly_derivative(Polynomial *p) {
    if (p->degree == 0) {
        Polynomial *d = create_polynomial(0);
        if (!d) return NULL;

        mpq_set_si(d->coeffs[0], 0, 1);
        return d;
    }

    Polynomial *d = create_polynomial(p->degree - 1);
    if (!d) return NULL;

    //(a_i * x_i)'=i * a_i * x^(i-1)
    mpq_t factor;
    mpq_init(factor);

    for (int i = 1; i <= p->degree; i++) {

        // factor = i (as rational i/1)
        mpq_set_si(factor, i, 1);

        // d[i-1] = factor * p->coeffs[i]
        mpq_mul(d->coeffs[i - 1], factor, p->coeffs[i]);
    }

    mpq_clear(factor);

    return d;
}

void poly_trim_degree(Polynomial *p) {

    // p[degree] = 0 -> degree--
    while (p->degree > 0) {
        if (mpz_sgn(mpq_numref(p->coeffs[p->degree])) != 0) break;
        p->degree--;
    }
}

/* Calculate rem(A,B) */
Polynomial* poly_remainder(Polynomial *A, Polynomial *B) {
    // B == 0
    if (!A || !B) return NULL;

    // division by zero polynomial
    if (B->degree == 0 && mpz_sgn(mpq_numref(B->coeffs[0])) == 0) {
        return NULL;
    }

    Polynomial *R = copy_polynomial(A);
    if (!R) return NULL;

    poly_trim_degree(R);

    if (R->degree < B->degree) {
        return R;
    }

    mpq_t c, tmp;
    mpq_init(c);
    mpq_init(tmp);

    while (R->degree >= B->degree) {

        // time start
        clock_t start = clock();

        int current_degree = R->degree;
        int k = R->degree - B->degree;

        // c = R_lc / B_lc 
        mpq_div(c, R->coeffs[R->degree], B->coeffs[B->degree]);

        // R = R - c * x^k * B 
        for (int i = 0; i <= B->degree; i++) {
            // tmp = c * B[i]
            mpq_mul(tmp, c, B->coeffs[i]);

            // R[i+k] -= tmp
            mpq_sub(R->coeffs[i + k], R->coeffs[i + k], tmp);
        }

        // // time end
        clock_t end = clock();

        // recalculate and update degree
        poly_trim_degree(R);


        double timeCost = (double)(end - start) / CLOCKS_PER_SEC;

        if (R->degree >= 0) {
            mpz_srcptr num = mpq_numref(R->coeffs[R->degree]);
            size_t bits = mpz_sizeinbase(num, 2);

            printf("Degree %d:\n", current_degree);
            printf("Time: %.6f s\n", timeCost);
            printf("Coeff numerator bits: %zu\n", bits);
        }

        printf("\n");
    }

    mpq_clear(c);
    mpq_clear(tmp);

    return R;
}

// (-1)*P
void poly_negative(Polynomial *p) {
    if (!p) return;

    int k = p->degree + 1;
    for (int i = 0; i < k; i++) {
        mpq_neg(p->coeffs[i], p->coeffs[i]);
    }
}

// Calculate P(value) with Hornor
// return sign of P(x): -1.0, 0.0, +1.0 
double poly_calculate_sign(Polynomial *p, mpq_t x) {
    if (!p) return 0.0;

    mpq_t acc;
    mpq_init(acc);
    mpq_set_si(acc, 0, 1);

    // method Hornor
    for (int i = p->degree; i >= 0; i--) {
        mpq_mul(acc, acc, x);
        mpq_add(acc, acc, p->coeffs[i]);
    }

    int s = mpz_sgn(mpq_numref(acc));

    mpq_clear(acc);

    if (s > 0) return 1.0;
    if (s < 0) return -1.0;
    return 0.0;
}

// Convert Polynomial to fmpq_poly_t
void polynomialConv(fmpq_poly_t poly, Polynomial *p)
{
    fmpq_poly_init(poly);

    for (int i = 0; i <= p->degree; i++) {
        fmpq_t c;
        fmpq_init(c);

        fmpq_set_mpq(c, p->coeffs[i]);
        fmpq_poly_set_coeff_fmpq(poly, i, c);

        fmpq_clear(c);
    }
}

// Calculate P(value) with FLINT
// return sign of P(x): -1.0, 0.0, +1.0 
int poly_sign_flint(fmpq_poly_t poly, mpq_t x)
{
    fmpq_t fx, res;

    fmpq_init(fx);
    fmpq_init(res);

    fmpq_set_mpq(fx, x);

    fmpq_poly_evaluate_fmpq(res, poly, fx);

    int s = fmpq_sgn(res);

    fmpq_clear(res);
    fmpq_clear(fx);

    return s;
}