#include "../include/pcca2.h"

int main(void) {
    srand((unsigned)time(NULL));

    // generate polynomial
    Polynomial*p=generate_random_polynomial(40,-100,100);
    
    printf("\n=== Test : random polynomial degree 40===\n");

    // sturm naif
    clock_t start1 = clock();
    mpq_t* b1 = sturm_naif(p,3);
    clock_t end1 = clock();
    
    double timeCost1 = (double)(end1 - start1) / CLOCKS_PER_SEC;
    printf("Time total cost for sturm_naif: %.6f s\n", timeCost1);

    if(verify_interval(b1,p)==0){
        printf("result of naif correct\n");
    }else{
        printf("result of naif wrong\n");
    }
    free_bound(b1);

    // sturm_habicht
    clock_t start = clock();
    mpq_t* b2 = sturm_habicht(p,3);
    clock_t end = clock();
    
    double timeCost = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time total cost for sturm_habicht: %.6f s\n", timeCost);

    if(verify_interval(b2,p)==0){
        printf("result of habicht correct\n");
    }else{
        printf("result of habicht wrong\n");
    }
    free_bound(b2);

    free_polynomial(p);



    return 0;
}