#include "../include/pcca2.h"

int main(void){
    Polynomial*p=generate_random_polynomial(40,-100,100);
    for(int i=2;i<=10;i++){
        clock_t start = clock();
        mpq_t* boundR = sturm_naif(p,i);
        clock_t end = clock();

        double timeCost = (double)(end - start) / CLOCKS_PER_SEC;
        printf("Number of divided intervals: %d ",i);
        printf("Time: %.6f s\n", timeCost);
        free_bound(boundR);
    }

    /*
    clock_t start = clock();
    mpq_t* boundR = sturm_naif(p,-1);
    clock_t end = clock();
    double timeCost = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Number of divided intervals = nbRacine ");
    printf("Time: %.6f s\n", timeCost);
    */

    free_polynomial(p);

}