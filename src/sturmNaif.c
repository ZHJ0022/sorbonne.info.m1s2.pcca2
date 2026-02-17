#include "../include/pcca2.h"

/*
    input: sturm suite(list of Polynomial) , number of the suite, upper and lower bounds of real roots
    output: bound[0] stores the number of real roots. For each real root i,
            its interval is given by bound[2*i - 1] and bound[2*i].
*/
mpq_t* bound_recu(Polynomial **sturmSuite, int nbSuite, mpq_t inf, mpq_t sup) {
    static int depth = 0;
    depth++;

    // allocate the 
    double *l1 = malloc(sizeof(double) * nbSuite);
    double *l2 = malloc(sizeof(double) * nbSuite);
    if (!l1 || !l2) {
        free(l1); free(l2);
        depth--;
        return NULL;
    }

    // Calculate the sign
    for (int i = 0; i < nbSuite; i++) {
        l1[i] = poly_calculate_sign(sturmSuite[i], inf);
        l2[i] = poly_calculate_sign(sturmSuite[i], sup);
    }

    // Calculate the number of sign changes
    int v1 = nb_sign_change(l1, nbSuite);
    int v2 = nb_sign_change(l2, nbSuite);
    int nbRacine = v1 - v2;

    free(l1);
    free(l2);

    if (nbRacine < 0) {
        printf("Error v1,v2: %d,%d\n", v1, v2);
        depth--;
        return NULL;
    }

    // allocate bound：size = 2*nbRacine + 1
    int size = 2 * nbRacine + 1;
    mpq_t *bound = malloc(sizeof(mpq_t) * size);
    if (!bound) {
        depth--;
        return NULL;
    }
    for (int i = 0; i < size; i++) {
        mpq_init(bound[i]);
        mpq_set_si(bound[i], 0, 1);
    }
    mpq_set_si(bound[0], nbRacine, 1);

    // root, return directly
    if (nbRacine == 0) {
        depth--;
        return bound;
    }

    // 1 root: Records the interval, return
    if (nbRacine == 1) {
        mpq_set(bound[1], inf);
        mpq_set(bound[2], sup);
        depth--;
        return bound;
    }

    // Multiple roots, recursively
    mpq_t diff, h, left, right, tmp, ii, n_q;
    mpq_init(diff);
    mpq_init(h);
    mpq_init(left);
    mpq_init(right);
    mpq_init(tmp);
    mpq_init(ii);
    mpq_init(n_q);

    int n = 2;
    mpq_sub(diff, sup, inf);        // diff = sup - inf
    mpq_set_si(n_q, n, 1);          // n_q=n=2
    mpq_div(h, diff, n_q);       // h = diff / n

    int index = 1;

    for (int i = 0; i < 2; i++) {
        // left = inf + i*h
        mpq_set(left, inf);
        mpq_set_si(ii, i, 1);
        mpq_mul(tmp, h, ii);
        mpq_add(left, left, tmp);

        // right = inf + (i+1)*h
        mpq_set(right, inf);
        mpq_set_si(ii, i + 1, 1);
        mpq_mul(tmp, h, ii);
        mpq_add(right, right, tmp);

        mpq_t *r = bound_recu(sturmSuite, nbSuite, left, right);
        if (r) {
            int rn = mpz_get_si(mpq_numref(r[0]));
            if (rn > 0) {
                for (int j = 0; j < rn; j++) {
                    if (index + 1 > nbRacine * 2) {
                        // printf("bound overflow error\n");
                        break;
                    }
                    mpq_set(bound[index],     r[1 + 2*j]);
                    mpq_set(bound[index + 1], r[2 + 2*j]);
                    index += 2;
                }
            }
            free_bound(r);
        }
    }

    mpq_clear(diff);
    mpq_clear(h);
    mpq_clear(left);
    mpq_clear(right);
    mpq_clear(tmp);
    mpq_clear(ii);
    mpq_clear(n_q);

    depth--;
    return bound;
}

/*
    input: polynomial p
    output: bound[0] stores the number of real roots. For each real root i,
            its interval is given by bound[2*i - 1] and bound[2*i].

*/ 

mpq_t* sturm_naif(Polynomial *p){
    int nbSuite = p->degree+1;
    // List of Polynomial objects representing the Sturm suite.
    Polynomial **sturmSuite = malloc(sizeof(Polynomial*) * nbSuite);
    if (!sturmSuite) return NULL;

    // Initialization the suite
    sturmSuite[0] = copy_polynomial(p);
    sturmSuite[1] = poly_derivative(p);   // p'

    if (!sturmSuite[0] || !sturmSuite[1]) {
        if (sturmSuite[0]) free_polynomial(sturmSuite[0]);
        if (sturmSuite[1]) free_polynomial(sturmSuite[1]);
        free(sturmSuite);
        return NULL;
    }

    // Calculate the suite
    for (int i = 1; i < p->degree; i++) {
        sturmSuite[i + 1] = poly_remainder(sturmSuite[i - 1], sturmSuite[i]);
        if (!sturmSuite[i + 1]) {
            for (int k = 0; k <= i; k++) free_polynomial(sturmSuite[k]);
            free(sturmSuite);
            return NULL;
        }

        // -rem()
        poly_negative(sturmSuite[i + 1]);

        // if the degree of the remainder = 0, stop calculate
        if (sturmSuite[i + 1]->degree == 0) {
            nbSuite = i + 2;
            break;
        }
    }

    // Cauchy bound: sup = bound, inf = -sup
    mpq_t sup, inf;
    mpq_init(sup);
    mpq_init(inf);

    cauchy_bound(sup, p);
    mpq_neg(inf, sup);
/*
    printf("cauchy_bound = ");
    gmp_printf("%Qd", inf);
    printf(", ");
    gmp_printf("%Qd", sup);
    printf("\n");
*/

    mpq_t *bound = bound_recu(sturmSuite, nbSuite, inf, sup);

    // free sturm suite
    for (int i = 0; i < nbSuite; i++) {
        free_polynomial(sturmSuite[i]);
    }
    free(sturmSuite);

    mpq_clear(sup);
    mpq_clear(inf);

    return bound;

}