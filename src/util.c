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

// Calculate the derivative of the polynomial
Polynomial* poly_derivative(Polynomial *p) {
    if (p->degree == 0) {
        Polynomial *d = create_polynomial(0);
        d->coeffs[0] = 0.0;
        return d;
    }

    Polynomial *d = create_polynomial(p->degree - 1);

    //(a_i * x_i)'=i * a_i * x^(i-1)
    for (int i = 1; i <= p->degree; i++) {
        d->coeffs[i - 1] = i * p->coeffs[i];
    }

    return d;
}

// Calculate rem(A,B)
Polynomial* poly_remainder(Polynomial *A,Polynomial *B) {
    Polynomial *R = copy_polynomial(A);

    if (A->degree < B->degree) {
        return R;
    }

    while (R->degree >= B->degree) {
        int k = R->degree - B->degree;
        double c = R->coeffs[R->degree] / B->coeffs[B->degree];

        // R = R - c * x^k * B
        for (int i = 0; i <= B->degree; i++) {
            R->coeffs[i + k] -= c * B->coeffs[i];
        }

        // update degree
        while (R->degree > 0 && R->coeffs[R->degree] == 0.0) {
            R->degree--;
        }
    }

    return R;
}

// Calculate P(value)
double poly_calculate(Polynomial *p, double value) {
    double result = 0.0;
    double x = 1.0;     //x_0

    for (int i = 0; i <= p->degree; i++) {
        result += p->coeffs[i] * x;
        x *= value;         //x_i+1
    }

    return result;
}