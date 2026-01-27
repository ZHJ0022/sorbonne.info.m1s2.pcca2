#include "../include/pcca2.h"

/*
    input: sturm suite(list of Polynomial) , number of the suite, upper and lower bounds of real roots
    output: bound[0] stores the number of real roots. For each real root i,
            its interval is given by bound[2*i - 1] and bound[2*i].
*/
double* bound_recu(Polynomial**sturmSuite, int nbSuite,double inf, double sup){

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

    int nbRacine=abs(v1-v2);


    double* bound=malloc(sizeof(double)*(nbRacine*2+1));
    bound[0]=nbRacine;

    free(l1);
    free(l2);

    // If there are no real roots, return directly.
    if(nbRacine==0){
        return bound;
    }

    /*  If there is only one real root
        -> complete the objective 
        -> record the current upper and lower bounds and then return.
     */ 
    if (nbRacine==1){
        bound[1]=inf;
        bound[2]=sup;
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
        sturmSuite[i+1]=poly_remainder(sturmSuite[i-1],sturmSuite[i]);

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
    
    double* bound=bound_recu(sturmSuite,nbSuite,inf,sup);
    for(int i=0;i<nbSuite;i++){
        free_polynomial(sturmSuite[i]);
    }
    free(sturmSuite);

    return bound;
}