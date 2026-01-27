#include "../include/pcca2.h"

int main(void) {
    // ---------- Test 1: (x+1)(x+2)(x+3) = x^3 + 6x^2 + 11x + 6 ----------
    Polynomial* p1 = create_polynomial(3);
    p1->coeffs[0] = 6;
    p1->coeffs[1] = 11;
    p1->coeffs[2] = 6;
    p1->coeffs[3] = 1;

    printf("=== Test 1: (x+1)(x+2)(x+3) ===\n");
    double* b1 = sturm_naif(p1);
    afficher_bound(b1);
    free(b1);
    free_polynomial(p1);

    // ---------- Test 2: (x-2)(x-4)(x-6)(x-8)(x-10) ----------
    /*
    Polynomial* p2 = create_polynomial(5);
    p2->coeffs[0] = -3840;
    p2->coeffs[1] = 4384;
    p2->coeffs[2] = -1800;
    p2->coeffs[3] = 340;
    p2->coeffs[4] = -30;
    p2->coeffs[5] = 1;

    printf("\n=== Test 2: (x-2)(x-4)(x-6)(x-8)(x-10) ===\n");
    */
    Polynomial* p2 = create_polynomial(2);
    p2->coeffs[0] = 1;
    p2->coeffs[1] = -2;
    p2->coeffs[2] = 1;
    printf("\n=== Test 2: (x-1)(x-1)===\n");
    double* b2 = sturm_naif(p2);
    afficher_bound(b2);
    free(b2);
    free_polynomial(p2);

    // ---------- Test 3: (x-1)(x-2)...(x-10) ----------
    Polynomial* p3 = create_polynomial(10);
    p3->coeffs[0]  = 3628800;
    p3->coeffs[1]  = -10628640;
    p3->coeffs[2]  = 12753576;
    p3->coeffs[3]  = -8409500;
    p3->coeffs[4]  = 3416930;
    p3->coeffs[5]  = -902055;
    p3->coeffs[6]  = 157773;
    p3->coeffs[7]  = -18150;
    p3->coeffs[8]  = 1320;
    p3->coeffs[9]  = -55;
    p3->coeffs[10] = 1;

    printf("\n=== Test 3: (x-1)(x-2)...(x-10) ===\n");
    double* b3 = sturm_naif(p3);
    afficher_bound(b3);
    free(b3);
    free_polynomial(p3);

    return 0;
}
