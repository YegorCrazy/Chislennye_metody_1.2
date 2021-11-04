#include <stdio.h>

typedef struct Matrix {
    long long m, n;
    long double **elem;
} matrix;

long double abs_d (long double x);

matrix *new_matrix (long long m, long long n);
matrix *copy_matrix (matrix *matr);
void delete_matrix (matrix *matr);
void print_matrix (matrix *matr);
void print_matrix_ext (matrix *A, matrix *f);
void read_matrix (FILE *input, matrix *matr);

void fill_matrix1 (matrix *A, matrix *f, long long n);
void fill_matrix2 (matrix *A, matrix *f, long long n);
long double q_M (long double M);
void fill_matrix3 (matrix *A, matrix *f, long long n);
void fill_matrix4 (matrix *A, matrix *f, long long n);
void fill_matrix5 (matrix *A, matrix *f, long long n);

void triangulate_matrix (matrix *matr, matrix *f);
long double determinant (matrix *matr);

matrix *reverse_matrix (matrix *A);
matrix *mult_matrix (matrix *A, matrix *B);
matrix *new_iter (matrix *x, matrix *A, matrix *f, long double omega);
long double diff (matrix *f, matrix *x);
