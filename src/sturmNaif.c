#include "../include/pcca2.h"

/*
    input: sturm suite(list of Polynomial) , number of the suite, upper and lower bounds of real roots
    output: bound[0] stores the number of real roots. For each real root i,
            its interval is given by bound[2*i - 1] and bound[2*i].
*/
double* bound_recu(Polynomial**sturmSuite, int nbSuite,double inf, double sup){
    static int depth = 0;
    depth++;
    /*
    if(depth > 2000){
        printf("depth overflow\n");
        return NULL;
    }
    */

    // printf("Enter [inf: %.15g, sup: %.15g], width: %.15g\n", inf, sup, sup - inf);

    double *l1=malloc(sizeof(double)*nbSuite);  // result of inf
    double *l2=malloc(sizeof(double)*nbSuite);  // result of sup

    // Calculate each value of l1, l2
    for(int i=0;i<nbSuite;i++){
        l1[i]=poly_calculate(sturmSuite[i],inf);
        l2[i]=poly_calculate(sturmSuite[i],sup);
    }
    
    // Calculate the number of sign changes for l1, l2
    int v1=nb_sign_change(l1,nbSuite);
    int v2=nb_sign_change(l2,nbSuite);

    int nbRacine=v1-v2;

    if(nbRacine<0){
        printf("Error v1,v2: %d,%d\n",v1,v2);
        
        printf("v1\n");
        v1 = nb_sign_changeDebug(l1,nbSuite);
        printf("v2\n");
        v1 = nb_sign_changeDebug(l2,nbSuite);
        
    }

    // printf("  -> v1: %d, v2: %d, nbRacine: %d\n", v1, v2, nbRacine);

    double* bound=malloc(sizeof(double)*(nbRacine*2+1));
    bound[0]=nbRacine;

    free(l1);
    free(l2);

    // If there are no real roots, return directly.
    if(nbRacine==0){
        depth--;
        return bound;
    }

    /*  If there is only one real root
        -> complete the objective 
        -> record the current upper and lower bounds and then return.
     */ 
    if (nbRacine==1){
        // printf("  SUCCESS: Isolated root in [%.15g, %.15g]\n", inf, sup);
        bound[1]=inf;
        bound[2]=sup;
        depth--;
        return bound;
    }

    /*

    If multiple real roots exist
    -> divide the range equally and continue the recursive operation.

    */
    double h=(sup-inf)/2;
    int index=1;                    // index of bound

    for(int i=0;i<2;i++){
        double* r=bound_recu(sturmSuite,nbSuite,inf+i*h,inf+(i+1)*h);
        if(r[0]>0.0){
            for(int j=0;j<(int)r[0];j++){
                if (index + 1 > nbRacine * 2) {
                    printf("bound overflow error\n");
                    printf("current interval[%.5f,%.5f], nb roots current: %f\n", inf+i*h,inf+(i+1)*h,r[0]);
                    printf("parent interval[%.5f,%.5f], nb roots parent: %d\n",inf,sup,nbRacine);
                    break;
                }
                bound[index]=r[1+j*2];
                bound[index+1]=r[(j+1)*2];
                index+=2;
            }
        }
        free(r);
    }

    return bound;
}

/*
    input: polynomial p
    output: bound[0] stores the number of real roots. For each real root i,
            its interval is given by bound[2*i - 1] and bound[2*i].

*/ 

double* sturm_naif(Polynomial* p){
    int nbSuite = p->degree+1;

    // List of Polynomial objects representing the Sturm suite.
    Polynomial** sturmSuite=malloc(sizeof(Polynomial*)*nbSuite);

    //Initialization the suite
    sturmSuite[0]=copy_polynomial(p);
    sturmSuite[1]=poly_derivative(p);   // p'

    //Calculate the suite
    for(int i=1;i<p->degree;i++){
        // Euclidean division

        // printf("Generating Sturm suite index %d...\n", i+1);
        sturmSuite[i+1]=poly_remainder(sturmSuite[i-1],sturmSuite[i]);

        // printf("->Remainder degree: %d\n", sturmSuite[i+1]->degree);

        // -rem()
        poly_negative(sturmSuite[i+1]);

        // if the degree of the remainder = 0, stop calculate
        if(sturmSuite[i+1]->degree==0){
            nbSuite=i+2;
            break;
        }
    }

    // cauchy bound
    double sup=cauchy_bound(p);
    double inf=(-1)*sup;
    printf("cauchy_bound = %.5f, %.5f\n", inf,sup);
    
    double* bound=bound_recu(sturmSuite,nbSuite,inf,sup);
    for(int i=0;i<nbSuite;i++){
        free_polynomial(sturmSuite[i]);
    }
    free(sturmSuite);

    return bound;
}