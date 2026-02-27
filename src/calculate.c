#include "../include/pcca2.h"

// Calculate cauchy bound
void cauchy_bound(mpq_t result, Polynomial *p) {
    // default result = 0
    mpq_set_si(result, 0, 1);

    poly_trim_degree(p);

    if (!p || p->degree == 0) {
        return;
    }

    int n = p->degree;

    mpq_t ratio, abs_ratio, maxv, one;
    mpq_init(ratio);
    mpq_init(abs_ratio);
    mpq_init(maxv);
    mpq_init(one);

    mpq_set_si(maxv, 0, 1);
    mpq_set_si(one, 1, 1);

    for (int i = 0; i < n; i++) {
        // ratio = a_i / a_n 
        mpq_div(ratio, p->coeffs[i], p->coeffs[n]);

        // abs_ratio = |ratio|
        mpq_set(abs_ratio, ratio);
        if (mpz_sgn(mpq_numref(abs_ratio)) < 0) {
            mpq_neg(abs_ratio, abs_ratio);
        }

        // maxv = max(maxv, abs_ratio)
        if (mpq_cmp(abs_ratio, maxv) > 0) {
            mpq_set(maxv, abs_ratio);
        }
    }

    // result = maxv + 1
    mpq_add(result, maxv, one);
    mpq_canonicalize(result);

    mpq_clear(ratio);
    mpq_clear(abs_ratio);
    mpq_clear(maxv);
    mpq_clear(one);
}

// Calculate the number of sign changes
// sign 1 positive, -1 negative, 0 zero
int nb_sign_change(double *l, int n) {
    int nb = 0;
    double curr,next;
    int index=0;
    while (index < n && l[index]==0.0) {
        index++;
    }
    if (index >= n) return 0;
    curr = l[index];

    for (int i = index+1; i < n; i++) {
        if (l[i] ==0.0){
            continue;
        }

        next=l[i];
        if(curr*next < 0.0){
            nb++;
        }
        curr=next;
    }

    return nb;
}

// Return 1 if the result is incorrect, and 0 if the result is correct.
int verify_interval(mpq_t *bound, Polynomial *p) {

    if (!bound || !p) return 1;

    // nb of roots
    int n = mpz_get_si(mpq_numref(bound[0]));

    for (int k = 0; k < p->nroots; k++) {

        int found = 0;

        for (int i = 1; i <= n; i++) {

            mpq_t *a = (mpq_t*)&bound[2*i - 1];
            mpq_t *b = (mpq_t*)&bound[2*i];

            // check if a <= root <= b 
            if (mpq_cmp(p->roots[k], *a) >= 0 &&
                mpq_cmp(p->roots[k], *b) <= 0) {

                found = 1;
                break;
            }
        }

        if (!found) return 1;
    }

    return 0;
}