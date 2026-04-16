#include "../include/pcca2.h"

void run_random_test(int case_id, int degree, long a, long b, int nbI) {
    printf("\n==================================================\n");
    printf("Random polynomial test #%d\n", case_id);
    printf("degree = %d, roots in [%ld, %ld], nbI = %d\n", degree, a, b, nbI);
    printf("==================================================\n");

    clock_t start_gen = clock();
    Polynomial *p = generate_random_polynomial(degree, a, b);
    clock_t end_gen = clock();

    double time_gen = (double)(end_gen - start_gen) / CLOCKS_PER_SEC;
    printf("Time to generate polynomial: %.6f s\n", time_gen);

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
    srand((unsigned)time(NULL));

    const int degree = 40;
    const long a = -100;
    const long b = 100;
    const int nbI = 3;
    const int ntests = 5;

    for (int i = 1; i <= ntests; i++) {
        run_random_test(i, degree, a, b, nbI);
    }

    return 0;
}