#include "../include/pcca2.h"

double* bound_recu(Polynomial**sturmSuite, int nbSuite,double inf, double sup){

    double *l1=malloc(sizeof(double)*nbSuite);  // result of inf
    double *l2=malloc(sizeof(double)*nbSuite);  // result of sup

    for(int i=0;i<nbSuite;i++){
        l1[i]=poly_calculate(sturmSuite[i],inf);
        l2[i]=poly_calculate(sturmSuite[i],sup);
    }
    
    int v1=nb_sign_change(l1,nbSuite);
    int v2=nb_sign_change(l2,nbSuite);

    int nbRacine=abs(v1-v2);

    double* bound=malloc(sizeof(double)*(nbRacine*2+1));
    bound[0]=nbRacine;

    if(nbRacine==0){
        return bound;
    }

    if (nbRacine==1){
        bound[1]=inf;
        bound[2]=sup;
        return bound;
    }

    double h=(sup-inf)/nbRacine;
    int index=1;                    // index of bound

    for(int i=0;i<nbRacine;i++){
        double* r=bound_recu(sturmSuite,nbSuite,inf+i*h,inf+(i+1)*h);
        if(r[0]>0){
            for(int j=0;j<r[0];j++){
                bound[index]=r[1+j*2];
                bound[index+1]=r[(j+1)*2];
                index+=2;
            }
        }
    }

    return bound;
}

double* sturm_naif(Polynomial* p){
    int nbSuite = p->degree+1;
    Polynomial** sturmSuite=malloc(sizeof(Polynomial*)*nbSuite);
    sturmSuite[0]=p;
    sturmSuite[1]=poly_derivative(p);

    for(int i=1;i<p->degree;i++){
        sturmSuite[i+1]=poly_remainder(sturmSuite[i-1],sturmSuite[i]);
        poly_negative(sturmSuite[i+1]);
    }

    // cauchy bound
    double sup=cauchy_bound(p);
    double inf=(-1)*sup;
    
    double* bound=bound_recu(sturmSuite,nbSuite,inf,sup);
    return bound;
}