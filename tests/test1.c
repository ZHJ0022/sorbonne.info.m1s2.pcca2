#include "../include/pcca2.h"

Polynomial* generate_structured_polynomial(int n) {
    if (n <= 0) return NULL;

    Polynomial *current = create_polynomial(0);
    if (!current) return NULL;
    mpq_set_si(current->coeffs[0], 1, 1);

    for (int k = 1; k <= n; k++) {
        // factor = x - k
        Polynomial *factor = create_polynomial(1);
        if (!factor) {
            free_polynomial(current);
            return NULL;
        }
        mpq_set_si(factor->coeffs[0], -k, 1);
        mpq_set_si(factor->coeffs[1], 1, 1);

        Polynomial *next = poly_mul_naive(current, factor);

        free_polynomial(current);
        free_polynomial(factor);

        if (!next) return NULL;
        current = next;
    }

    // store roots for verify_interval
    current->nroots = n;
    current->roots = (mpq_t*)malloc(sizeof(mpq_t) * n);
    if (!current->roots) {
        free_polynomial(current);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        mpq_init(current->roots[i]);
        mpq_set_si(current->roots[i], i + 1, 1);
    }

    return current;
}

void run_structured_test(int n, int nbI) {
    printf("\n==================================================\n");
    printf("Structured polynomial: P_%d(x) = (x-1)(x-2)...(x-%d)\n", n, n);
    printf("degree = %d\n", n);
    printf("==================================================\n");

    Polynomial *p = generate_structured_polynomial(n);

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
    int nbI = 3;

    int degrees[] = {5, 10, 15, 20, 25, 30};
    int ncases = sizeof(degrees) / sizeof(degrees[0]);

    for (int i = 0; i < ncases; i++) {
        run_structured_test(degrees[i], nbI);
    }

    return 0;
}