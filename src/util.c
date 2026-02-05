#include "../include/pcca2.h"

// Create a polynomial
Polynomial* create_polynomial(int degree) {
    Polynomial *p = (Polynomial*)malloc(sizeof(Polynomial));
    if (!p) return NULL;

    p->degree = degree;
    p->coeffs = (double*)calloc(degree + 1, sizeof(double));
    if (!p->coeffs) {
        free(p);
        return NULL;
    }

    p->roots = NULL;
    p->nroots = 0;

    return p;
}

// Free polynomial memory
void free_polynomial(Polynomial *p) {
    if (p) {
        if (p->coeffs) free(p->coeffs);
        if (p->roots)  free(p->roots);
        free(p);
    }
}

// Copy a polynomial
Polynomial* copy_polynomial(const Polynomial *p) {
    if (!p) return NULL;

    Polynomial *copy = create_polynomial(p->degree);
    if (!copy) return NULL;

    memcpy(copy->coeffs, p->coeffs, (p->degree + 1) * sizeof(double));

    copy->nroots = p->nroots;
    if (p->nroots > 0 && p->roots) {
        copy->roots = (double*)malloc(p->nroots * sizeof(double));

        if (!copy->roots) {
            free_polynomial(copy);
            return NULL;
        }

        memcpy(copy->roots, p->roots, p->nroots * sizeof(double));
    }

    return copy;
}

Polynomial* poly_mul_naive(Polynomial *A, Polynomial *B) {
    Polynomial *C = create_polynomial(A->degree + B->degree);
    if (!C) return NULL;

    for (int i = 0; i <= A->degree; i++) {
        for (int j = 0; j <= B->degree; j++) {
            C->coeffs[i + j] += A->coeffs[i] * B->coeffs[j];
        }
    }
    return C;
}

// Generate polynomial from n random real roots in [a,b]:
// P(x) = (x - r_1)(x - r_2)...(x - r_n), degree = n
Polynomial* generate_random_polynomial(int n, double a, double b) {
    if (n < 0) return NULL;

    Polynomial *p = create_polynomial(n);
    if (!p) return NULL;

    // allocate roots
    p->roots = (double*)malloc(n * sizeof(double));
    if (n > 0 && !p->roots) {
        free_polynomial(p);
        return NULL;
    }
    p->nroots = n;

    // temp stores the calculation results and updates p.
    double *temp = (double*)malloc((n + 1) * sizeof(double));
    if (!temp) {
        free_polynomial(p);
        return NULL;
    }

    p->coeffs[0] = 1.0;

    for (int k = 0; k < n; k++) {
        double r = a + ((double)rand() / RAND_MAX) * (b - a);
        p->roots[k] = r;

        // Initialize temp(result) to 0
        for (int i = 0; i <= n; i++) {
            temp[i] = 0.0;
        }

        for (int i = 0; i <= k; i++) {
            temp[i] += (-r) * p->coeffs[i];
        }

        for (int i = 0; i <= k; i++) {
            temp[i + 1] += 1 * p->coeffs[i];
        }

        // copy the result
        for (int i = 0; i <= k + 1; i++) {
            p->coeffs[i] = temp[i];
        }
    }

    free(temp);

    return p;
}


// afficher the result of the sturm
void afficher_bound(const double* bound) {
    if (!bound) {
        printf("bound = NULL\n");
        return;
    }

    int n = (int)bound[0];
    printf("Nombre de racines reelles (selon Sturm) = %d\n", n);

    for (int i = 1; i <= n; i++) {
        double a = bound[2*i - 1];
        double b = bound[2*i];
        printf("  Racine %d : [%.15g, %.15g]\n", i, a, b);
    }
}

// Calculate the derivative of the polynomial
Polynomial* poly_derivative(Polynomial *p) {
    if (p->degree == 0) {
        Polynomial *d = create_polynomial(0);
        d->coeffs[0] = 0.0;
        return d;
    }

    Polynomial *d = create_polynomial(p->degree - 1);

    //(a_i * x_i)'=i * a_i * x^(i-1)
    for (int i = 1; i <= p->degree; i++) {
        d->coeffs[i - 1] = i * p->coeffs[i];
    }

    return d;
}

// Calculate rem(A,B)
Polynomial* poly_remainder(Polynomial *A,Polynomial *B) {
    if (B->degree == 0) {
        return create_polynomial(0);
    }

    Polynomial *R = copy_polynomial(A);

    if (A->degree < B->degree) {
        return R;
    }

    while (R->degree >= B->degree) {
        int k = R->degree - B->degree;
        double c = R->coeffs[R->degree] / B->coeffs[B->degree];

        // R = R - c * x^k * B
        for (int i = 0; i <= B->degree; i++) {
            if (i + k > A->degree){
                printf("remainder overflow error");
            }
            R->coeffs[i + k] -= c * B->coeffs[i];
        }

        // update degree
        while (R->degree > 0 && fabs(R->coeffs[R->degree]) < EPSILON) {
            R->degree--;
        }
    }

    return R;
}

void poly_negative(Polynomial *p){
    int k=p->degree+1;
    for(int i=0;i<k;i++){
        p->coeffs[i]*=-1;
    }
}

// Calculate P(value)
double poly_calculate(Polynomial *p, double value) {
    double result = 0.0;
    double x = 1.0;     //x_0

    for (int i = 0; i <= p->degree; i++) {
        result += p->coeffs[i] * x;
        x *= value;         //x_i+1
    }

    return result;
}

// Calculate cauchy bound
double cauchy_bound(Polynomial*p){
    if (p->degree == 0) {
        return 0.0;
    }

    double max=0.0;
    int n=p->degree;

    for (int i = 0; i < n; i++) {
        double r = fabs(p->coeffs[i] / p->coeffs[n]);
        if (r > max) {
            max=r;
        }
    }

    return max+1;
}

// Calculate the number of sign changes
int nb_sign_change(double *l, int n) {
    int nb = 0;

    double curr,next;
    int index=0;
    curr=l[index];

    while (index < n && l[index] ==0.0 ) {
        index++;
    }
    if (index >= n) return 0;

    curr = l[index];

    for (int i = index+1; i < n; i++) {
        if (l[i] == 0.0){
            continue;
        }

        next=l[i];
        if(curr*next<0){
            nb++;
        }
        curr=next;

    }

    return nb;
}

// Return 1 if the result is incorrect, and 0 if the result is correct.
int verify_interval(double *bound, Polynomial *p){
    int n = (int)bound[0];

    for (int k = 0; k < p->nroots; k++) {
        double r = p->roots[k];
        int found = 0;

        for (int i = 1; i <= n; i++) {
            double a = bound[2*i - 1];
            double b = bound[2*i];
            if (r >= a && r <= b){
                found = 1; break;
            }
        }
        if (!found) return 1;
    }
    return 0;
}