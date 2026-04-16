#include "../include/pcca2.h"

void set_mpq_pow2(mpq_t q, unsigned int k) {
    mpz_t z;
    mpz_init(z);
    mpz_ui_pow_ui(z, 2, k);  
    mpq_set_z(q, z);
    mpq_canonicalize(q);
    mpz_clear(z);
}

void set_mpq_root(mpq_t q, unsigned long num, unsigned int k) {
    mpz_t den, n;
    mpz_init(den);
    mpz_init_set_ui(n, num);
    mpz_ui_pow_ui(den, 2, k);

    mpq_set_num(q, n);
    mpq_set_den(q, den);
    mpq_canonicalize(q);

    mpz_clear(den);
    mpz_clear(n);
}

Polynomial* generate_close_roots_integer_polynomial(int degree, unsigned int k) {
    if (degree <= 0) return NULL;

    Polynomial *current = create_polynomial(0);
    if (!current) return NULL;
    mpq_set_si(current->coeffs[0], 1, 1);

    for (int i = 0; i < degree; i++) {
        Polynomial *factor = create_polynomial(1);
        if (!factor) {
            free_polynomial(current);
            return NULL;
        }

        
        // coeffs[0] = -i 
        // mpq_set_si(factor->coeffs[0], -i, 1);

        // coeffs[1] = 2^k 
        // set_mpq_pow2(factor->coeffs[1], k);
        
        mpq_neg(factor->coeffs[0], factor->coeffs[0]);   // coeffs[0] = -i / 2^k
        mpq_set_si(factor->coeffs[1], 1, 1);             // coeffs[1] = 1

        Polynomial *next = poly_mul_naive(current, factor);

        free_polynomial(current);
        free_polynomial(factor);

        if (!next) return NULL;
        current = next;
    }

    // store roots for verify_interval
    current->nroots = degree;
    current->roots = (mpq_t*)malloc(sizeof(mpq_t) * degree);
    if (!current->roots) {
        free_polynomial(current);
        return NULL;
    }

    for (int i = 0; i < degree; i++) {
        mpq_init(current->roots[i]);
        set_mpq_root(current->roots[i], (unsigned long)i, k);
    }

    return current;
}

void run_close_root_test(int case_id, int degree, unsigned int k, int nbI) {
    printf("\n==================================================\n");
    printf("Close roots test #%d\n", case_id);
    printf("Polynomial: prod_{i=0}^{%d} (2^%u x - i)\n", degree - 1, k);
    printf("Degree = %d, roots = { i / 2^%u | 0 <= i <= %d }, nbI = %d\n",
           degree, k, degree - 1, nbI);
    printf("==================================================\n");

    Polynomial *p = generate_close_roots_integer_polynomial(degree, k);


    // sturm naif
    printf("\nNaive Sturm\n");
    mpq_t *b1 = sturm_naif(p, nbI);
    if (verify_interval(b1, p) == 0)
        printf("Naive Sturm: result correct\n");
    else
        printf("Naive Sturm: result wrong\n");


    free_bound(b1);

    //sturm habicht
    printf("\nSturm-Habicht\n");
    mpq_t *b2 = sturm_habicht(p, nbI);
    if (verify_interval(b2, p) == 0)
        printf("Sturm-Habicht: result correct\n");
    else
        printf("Sturm-Habicht: result wrong\n");

    free_bound(b2);

    free_polynomial(p);
}

int main(void) {

    const int degree = 40;
    const int nbI = 3;

    run_close_root_test(1, degree, 16, nbI);
    run_close_root_test(2, degree, 24, nbI);
    run_close_root_test(3, degree, 32, nbI);

    return 0;
}