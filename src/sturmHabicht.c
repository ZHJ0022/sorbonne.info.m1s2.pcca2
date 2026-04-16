#include "../include/pcca2.h"

Polynomial** create_sturm_habicht(Polynomial *P)
{
    if (P == NULL)
        return NULL;

    if (poly_is_zero(P))
        return NULL;

    int p = P->degree;

    Polynomial **seq = (Polynomial**)calloc(p + 1, sizeof(Polynomial*));
    if (seq == NULL)
        return NULL;

    mpq_t lcP;
    mpq_init(lcP);
    poly_get_lc(lcP, P);

    // Q = 1, P'Q = P'
    Polynomial *dP = poly_derivative(P);

    // case 1 : lcof(P) = 1
    if (poly_is_monic(P)) {
        Polynomial *R = poly_remainder(dP, P);
        poly_trim_degree(R);

        int r = poly_is_zero(R) ? -1 : R->degree;

        // StHa_p(P) = P
        seq[p] = copy_polynomial(P);

        // StHa_{p-1}(P) = R
        if (p - 1 >= 0)
            seq[p - 1] = copy_polynomial(R);

        if (r >= 0) {
            // StHa_r(P) = lcof(R)^(p-r+1) * R
            mpq_t lcR, pow_lcR;
            mpq_init(lcR);
            mpq_init(pow_lcR);

            poly_get_lc(lcR, R);

            mpz_t num, den;
            mpz_init(num);
            mpz_init(den);

            mpz_pow_ui(num, mpq_numref(lcR), (unsigned long)(p - r + 1));
            mpz_pow_ui(den, mpq_denref(lcR), (unsigned long)(p - r + 1));
            mpq_set_num(pow_lcR, num);
            mpq_set_den(pow_lcR, den);
            mpq_canonicalize(pow_lcR);

            mpz_clear(num);
            mpz_clear(den);

            if (seq[r] != NULL)
                free_polynomial(seq[r]);

            seq[r] = copy_polynomial(R);
            poly_scalar_mul(seq[r], pow_lcR);

            mpq_clear(lcR);
            mpq_clear(pow_lcR);

            // StHa_j(P) = 0, for r+1 <= j <= p-2
            for (int j = r + 1; j <= p - 2; j++) {
                if (seq[j] == NULL) {
                    seq[j] = create_polynomial(0);
                    mpq_set_ui(seq[j]->coeffs[0], 0, 1);
                    seq[j]->degree = 0;
                }
            }

            // StHa_j(P) = δ_{p-j-1} * Sres_j(P,p,R,r)
            for (int j = 0; j < r; j++) {
                seq[j] = subresultant_j(P, p, R, r, j);

                if (delta_sign(p - j - 1) < 0)
                    poly_negative(seq[j]);

                poly_trim_degree(seq[j]);
            }
        } else {
            
            // R=0
            for (int j = 0; j <= p - 2; j++) {
                if (seq[j] == NULL) {
                    seq[j] = create_polynomial(0);
                    mpq_set_ui(seq[j]->coeffs[0], 0, 1);
                    seq[j]->degree = 0;
                }
            }
        }

        free_polynomial(R);
        free_polynomial(dP);
        mpq_clear(lcP);

        return seq;
    }

    // case 2 : lcof(P) != 1
    // Q = 1,q = 0

    // StHa_p(P) = lcof(P) * P
    seq[p] = copy_polynomial(P);
    poly_scalar_mul(seq[p], lcP);

    // StHa_{p-1}(P) = lcof(P) * P'
    if (p - 1 >= 0) {
        seq[p - 1] = copy_polynomial(dP);
        poly_scalar_mul(seq[p - 1], lcP);
    }

    // StHa_j(P) = δ_{p-j-1} * Sres_j(P,p,P',p-1) / lcof(P), j < p-1
    for (int j = 0; j < p - 1; j++) {
        seq[j] = subresultant_j(P, p, dP, p - 1, j);

        if (delta_sign(p - j - 1) < 0)
            poly_negative(seq[j]);

        if (poly_scalar_div(seq[j], lcP) != 0) {
            for (int k = 0; k <= p; k++) {
                if (seq[k] != NULL)
                    free_polynomial(seq[k]);
            }
            free(seq);

            free_polynomial(dP);
            mpq_clear(lcP);
            return NULL;
        }

        poly_trim_degree(seq[j]);
    }

    free_polynomial(dP);
    mpq_clear(lcP);

    return seq;
}

mpq_t* sturm_habicht(Polynomial *p, int nbI){
    int nbSuite = p->degree +1;

    clock_t start = clock();
    Polynomial **sturmHabicht= create_sturm_habicht(p);
    clock_t end = clock();
    
    double timeCost = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time cost for calculate sequence sturm_habicht: %.6f s\n", timeCost);

    if (!sturmHabicht) return NULL;

    // Cauchy bound: sup = bound, inf = -sup
    mpq_t sup, inf;
    mpq_init(sup);
    mpq_init(inf);

    cauchy_bound(sup, p);
    mpq_neg(inf, sup);

    clock_t start2 = clock();

    mpq_t *bound = bound_recu(sturmHabicht, nbSuite, inf, sup, nbI);
    clock_t end2 = clock();

    double timeCost2 = (double)(end2 - start2) / CLOCKS_PER_SEC;
    printf("Time cost for calculate interval sturm_habicht: %.6f s\n", timeCost2);


    // free sturm_habicht
    for (int i = 0; i < nbSuite; i++) {
        if (sturmHabicht[i] != NULL)
            free_polynomial(sturmHabicht[i]);
    }
    free(sturmHabicht);

    mpq_clear(sup);
    mpq_clear(inf);

    return bound;
}