#include "../include/pcca2.h"

// create a new matrix of mpq_t
mpq_t** alloc_matrix(int rows, int cols)
{
    mpq_t **M = malloc(rows * sizeof(mpq_t*));
    if (!M) return NULL;

    for (int i = 0; i < rows; i++) {
        M[i] = malloc(cols * sizeof(mpq_t));
        if (!M[i]) {
            for (int k = 0; k < i; k++) {
                for (int j = 0; j < cols; j++) mpq_clear(M[k][j]);
                free(M[k]);
            }
            free(M);
            return NULL;
        }
        for (int j = 0; j < cols; j++) {
            mpq_init(M[i][j]);
            mpq_set_ui(M[i][j], 0, 1);
        }
    }
    return M;
}

// free the matrix of mpq_t
void free_matrix(mpq_t **M, int rows, int cols)
{
    if (!M) return;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) mpq_clear(M[i][j]);
        free(M[i]);
    }
    free(M);
}

// create Sylvester matrix of P and S
mpq_t** create_sylv_j_matrix(Polynomial *P, Polynomial *S, int j,
                            int *out_rows, int *out_cols)
{
    int p = P->degree;
    int s = S->degree;

    int rowsP = s - j;
    int rowsS = p - j;
    int rows = rowsP + rowsS;
    int cols = p + s - j;

    if (j < 0 || j >= (p < s ? p : s)) return NULL;

    mpq_t **M = alloc_matrix(rows, cols);
    if (!M) return NULL;

    // P block: X^(s-j-1)P, ..., XP, P
    for (int r = 0; r < rowsP; r++) {
        int shift = rowsP - 1 - r;   // from high shift down to 0 

        for (int k = 0; k <= p; k++) {
            int deg = k + shift;                     // actual monomial degree 
            int col = (p + s - j - 1) - deg;        // basis is descending 
            mpq_set(M[r][col], P->coeffs[k]);
        }
    }

    // S block: X^(p-j-1)S, ..., XS, S
    for (int r = 0; r < rowsS; r++) {
        int shift = rowsS - 1 - r;
        int row = rowsP + r;

        for (int k = 0; k <= s; k++) {
            int deg = k + shift;
            int col = (p + s - j - 1) - deg;
            mpq_set(M[row][col], S->coeffs[k]);
        }
    }

    if (out_rows) *out_rows = rows;
    if (out_cols) *out_cols = cols;
    return M;
}

// convert mpq_t(GMP) to fmpq_t(FLINT)
static void mpq_to_fmpq(fmpq_t y, const mpq_t x)
{
    fmpz_set_mpz(fmpq_numref(y), mpq_numref(x));
    fmpz_set_mpz(fmpq_denref(y), mpq_denref(x));
    fmpq_canonicalise(y);
}

// convert fmpq_t(FLINT) to mpq_t(GMP)
static void fmpq_to_mpq(mpq_t y, const fmpq_t x)
{
    mpz_set(mpq_numref(y), fmpz_promote_val(fmpq_numref(x)));
    mpz_set(mpq_denref(y), fmpz_promote_val(fmpq_denref(x)));
    mpq_canonicalize(y);
}

void sylv_j_det(mpq_t out, mpq_t **M, int rows, int cols, int d)
{
    fmpq_mat_t A;
    fmpq_t det;
    fmpq_init(det);
    fmpq_mat_init(A, rows, rows);

    for (int i = 0; i < rows; i++) {
        for (int c = 0; c < rows - 1; c++) {
            mpq_to_fmpq(fmpq_mat_entry(A, i, c), M[i][c]);
        }
        mpq_to_fmpq(fmpq_mat_entry(A, i, rows - 1), M[i][cols - d - 1]);
    }

    fmpq_mat_det(det, A);
    fmpq_to_mpq(out, det);

    fmpq_clear(det);
    fmpq_mat_clear(A);
}

Polynomial* subresultant_j(Polynomial *P, Polynomial *S, int j)
{
    if (!P || !S) return NULL;

    int p = P->degree;
    int s = S->degree;
    int minps = (p < s ? p : s);
    if (j < 0 || j >= minps) return NULL;

    int rows, cols;
    mpq_t **M = create_sylv_j_matrix(P, S, j, &rows, &cols);
    if (!M) return NULL;

    Polynomial *R = create_polynomial(j);
    if (!R) {
        free_matrix(M, rows, cols);
        return NULL;
    }

    for (int d = 0; d <= j; d++) {
        sylv_j_coefficient_det(R->coeffs[d], M, rows, cols, d);
    }

    poly_trim_degree(R);
    free_matrix(M, rows, cols);
    return R;
}