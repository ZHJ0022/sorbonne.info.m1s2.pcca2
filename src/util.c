#include "../include/pcca2.h"

// Create a polynomial
Polynomial* create_polynomial(int degree) {
    Polynomial *p = (Polynomial*)malloc(sizeof(Polynomial));
    if (!p) return NULL;

    p->degree = degree;

    // allocate coeff
    p->coeffs = (mpq_t*)malloc((degree + 1) * sizeof(mpq_t));
    if (!p->coeffs) {
        free(p);
        return NULL;
    }

    // init coeff
    for (int i = 0; i <= degree; i++) {
        mpq_init(p->coeffs[i]);
        mpq_set_si(p->coeffs[i], 0, 1);   /* 0/1 */
    }

    p->roots = NULL;
    p->nroots = 0;

    return p;
}

// Free polynomial memory
void free_polynomial(Polynomial *p) {
    if (!p) return;

    if (p->coeffs) {
        for (int i = 0; i <= p->degree; i++) {
            mpq_clear(p->coeffs[i]);
        }
        free(p->coeffs);
    }

    if (p->roots) {
        for (int i = 0; i < p->nroots; i++) {
            mpq_clear(p->roots[i]);
        }
        free(p->roots);
    }

    free(p);
}

// Free mpq_t*
void free_bound(mpq_t *bound) {
    if (!bound) return;
    int n = mpz_get_si(mpq_numref(bound[0]));
    int size = 2*n + 1;
    for (int i = 0; i < size; i++) mpq_clear(bound[i]);
    free(bound);
}

// Copy a polynomial
Polynomial* copy_polynomial(Polynomial *p) {
    if (!p) return NULL;

    Polynomial *copy = create_polynomial(p->degree);
    if (!copy) return NULL;

    //copy coeff
    for (int i = 0; i <= p->degree; i++) {
        mpq_set(copy->coeffs[i], p->coeffs[i]);
    }

    //copy root
    copy->nroots = p->nroots;
    copy->roots = NULL;

    if (p->nroots > 0 && p->roots) {
        copy->roots = (mpq_t*)malloc(p->nroots * sizeof(mpq_t));
        if (!copy->roots) {
            free_polynomial(copy);
            return NULL;
        }

        for (int i = 0; i < p->nroots; i++) {
            mpq_init(copy->roots[i]);
            mpq_set(copy->roots[i], p->roots[i]);
        }
    }

    return copy;
}

// Generate polynomial from n random real roots in [a,b]:
// P(x) = (x - r_1)(x - r_2)...(x - r_n), degree = n
Polynomial* generate_random_polynomial(int n, long a, long b) {
    if (n < 0) return NULL;
    if (n > 0 && a > b) return NULL;

    Polynomial *p = create_polynomial(n);
    if (!p) return NULL;

    p->nroots = n;
    p->roots = NULL;
    if (n > 0) {
        p->roots = (mpq_t*)malloc(n * sizeof(mpq_t));
        if (!p->roots) {
            free_polynomial(p);
            return NULL;
        }
        for (int i = 0; i < n; i++) {
            mpq_init(p->roots[i]);
            mpq_set_si(p->roots[i], 0, 1);
        }
    }

    // temp polynomial coefficients (mpq) size n+1
    mpq_t *temp = (mpq_t*)malloc((n + 1) * sizeof(mpq_t));
    if (!temp) {
        free_polynomial(p);
        return NULL;
    }
    for (int i = 0; i <= n; i++) {
        mpq_init(temp[i]);
        mpq_set_si(temp[i], 0, 1);
    }

    // init: p(x) = 1
    mpq_set_si(p->coeffs[0], 1, 1);

    mpq_t rq, neg_rq, prod;
    mpq_init(rq);
    mpq_init(neg_rq);
    mpq_init(prod);
    long max_denominator = 10;

    for (int k = 0; k < n; k++) {

        // Randomly generate r (numerator) and q (denominator)
        long q = 1 + rand() % max_denominator;
        long r = a*q + rand() % ((b-a)*q + 1);
        mpq_set_si(rq, r, q);
        mpq_canonicalize(rq);

        mpq_set(p->roots[k], rq);

        /* neg_rq = -r */
        mpq_neg(neg_rq, rq);

        // temp[...] = 0
        for (int i = 0; i <= n; i++) {
            mpq_set_si(temp[i], 0, 1);
        }

        // Multiply current poly(degree k) by (x - r):
        // new[i]   += (-r) * old[i]
        // new[i+1] += 1    * old[i]
        for (int i = 0; i <= k; i++) {
            /* prod = (-r) * p->coeffs[i] */
            mpq_mul(prod, neg_rq, p->coeffs[i]);
            mpq_add(temp[i], temp[i], prod);

            // temp[i+1] += p->coeffs[i]
            mpq_add(temp[i + 1], temp[i + 1], p->coeffs[i]);
        }

        // copy back to p coeffs up to degree k+1
        for (int i = 0; i <= k + 1; i++) {
            mpq_set(p->coeffs[i], temp[i]);
        }
    }

    /* cleanup temp */
    mpq_clear(rq);
    mpq_clear(neg_rq);
    mpq_clear(prod);

    for(int i = 0; i <= n; i++) {
        mpq_clear(temp[i]);
    }
    free(temp);

    // gmp_printf("Coeff degree n: %Qd\n", p->coeffs[n]);

    return p;
}

// afficher the result of the sturm
void afficher_bound(mpq_t *bound) {
    if (!bound) {
        printf("bound = NULL\n");
        return;
    }

    // bound[0] nb roots
    int n = mpz_get_si(mpq_numref(bound[0]));

    printf("Nombre de racines reelles (selon Sturm) = %d\n", n);

    for (int i = 1; i <= n; i++) {
        printf("  Racine %d : [", i);
        gmp_printf("%Qd", bound[2*i - 1]);
        printf(", ");
        gmp_printf("%Qd", bound[2*i]);
        printf("]\n");
    }
}
