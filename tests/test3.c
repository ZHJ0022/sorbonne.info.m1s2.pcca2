#include "../include/pcca2.h"

int main(void){

        clock_t start1 = clock();
        Polynomial*p=generate_random_polynomial(40,-100,100);
        clock_t end1 = clock();

        double time1 = (double)(end1 - start1) / CLOCKS_PER_SEC;
        printf("Time generate polynomial: %.6f s\n",time1);

        int nbSuite = p->degree+1;

        Polynomial **sturmSuite = malloc(sizeof(Polynomial*) * nbSuite);

        sturmSuite[0] = copy_polynomial(p);
        sturmSuite[1] = poly_derivative(p);  

        for (int i = 1; i < p->degree; i++) {

            printf("Sturm suite [%d]: \n",i);
            sturmSuite[i + 1] = poly_remainder(sturmSuite[i - 1], sturmSuite[i]);

            // -rem()
            poly_negative(sturmSuite[i + 1]);

            // if the degree of the remainder = 0, stop calculate
            if (sturmSuite[i + 1]->degree == 0) {
                nbSuite = i + 2;
                break;
            }
        }

        // free sturm suite
        for (int i = 0; i < nbSuite; i++) {
            free_polynomial(sturmSuite[i]);
        }
        free(sturmSuite);
        free_polynomial(p);


}