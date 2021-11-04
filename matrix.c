#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct Matrix {
    long long m, n; //матрица, m строк, n столбцов
    long double **elem;
} matrix;

long double abs_d (long double x) {
    if (x < 0) return -x;
    return x;
}

matrix *new_matrix (long long m, long long n) { //создание матрицы
    matrix *res = malloc(sizeof(matrix));
    if (res == NULL) return res;
    res->m = m;
    res->n = n;
    long double **elem = malloc(sizeof(long double *) * n);
    for (long long i = 0; i < n; ++i) {
        elem[i] = calloc(m, sizeof(long double)); // память по столбцам для удобного выбора ведущего элемента
    }
    res->elem = elem;
    return res;
}

matrix *copy_matrix (matrix *matr) {
    long long m = matr->m;
    long long n = matr->n;
    matrix *new = new_matrix(matr->m, matr->n);
    for (long long i = 0; i < n; ++i) {
        for (long long j = 0; j < m; ++j) {
            (new->elem)[i][j] = (matr->elem)[i][j];
        }
    }
    return new;
}

void delete_matrix (matrix *matr) { //удаление матрицы
    long long n = matr->n;
    for (long long i = 0; i < n; ++i) {
        free((matr->elem)[i]);
    }
    free(matr->elem);
    free(matr);
}

void print_matrix (matrix *matr) { //вывод матрицы
    long long m = matr->m;
    long long n = matr->n;
    for (long long i = 0; i < m; ++i) {
        for (long long j = 0; j < n; ++j) {
            printf("%Lf\t", (matr->elem)[j][i]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_matrix_ext (matrix *A, matrix *f) { //вывод матрицы
    long long m = A->m;
    long long n = A->n;
    for (long long i = 0; i < m; ++i) {
        for (long long j = 0; j < n; ++j) {
            printf("%Lf\t", (A->elem)[j][i]);
        }
        printf("%Lf\n", (f->elem)[0][i]);
    }
    printf("\n");
}

void read_matrix (FILE *input, matrix *matr) { //чтение матрицы из файла
    long long m = matr->m;    
    long long n = matr->n;
    for (long long i = 0; i < m; ++i) {
        for (long long j = 0; j < n; ++j) {
            fscanf(input, "%Lf", &((matr->elem)[j][i]));
        }
    }
}

void fill_matrix1 (matrix *A, matrix *f, long long n) {
    if (n != 20) {
        printf("ERROR ERROR ERROR\n");
        return;
    }
    long double m1 = 8, n1 = 20;
    for (long long i = 0; i < n; ++i) {
        for (long long j = 0; j < n; ++j) {
            if (i != j) {
                (A->elem)[i][j] = (i + j + 2) / (m1 + n1);
            } else {
                (A->elem)[i][j] = n1 + (m1 * m1) + ((j + 1) / m1) + ((i + 1) / n1);
            }
        }
    }
    for (long long i = 0; i < n; ++i) {
        (f->elem)[0][i] = 200 + (50 * (i + 1));
    }
}

void fill_matrix2 (matrix *A, matrix *f, long long n) {
    if (n != 25) {
        printf("ERROR ERROR ERROR\n");
        return;
    }
    long double m1 = 10, n1 = 25;
    for (long long i = 0; i < n; ++i) {
        for (long long j = 0; j < n; ++j) {
            if (i != j) {
                (A->elem)[i][j] = (i + j + 2) / (m1 + n1);
            } else {
                (A->elem)[i][j] = n1 + (m1 * m1) + ((j + 1) / m1) + ((i + 1) / n1);
            }
        }
    }
    for (long long i = 0; i < n; ++i) {
        (f->elem)[0][i] = (i + 1) * (i + 1) - n1;
    }
}

long double q_M (long double M) {
    return 1.001 - (2 * M * 0.001);
}

void fill_matrix3 (matrix *A, matrix *f, long long n) {
    if (n != 100) {
        printf("ERROR ERROR ERROR\n");
        return;
    }
    long double M1 = 4, n1 = 100;
    long double qM = q_M(M1);
    for (long long i = 0; i < n; ++i) {
        for (long long j = 0; j < n; ++j) {
            if (i != j) {
                (A->elem)[i][j] = pow(qM, (double)(i + j + 2)) + 0.1 * ((double)j - (double)i);
            } else {
                (A->elem)[i][j] = pow(abs_d(qM - 1), (double)(i + j + 2));
                if ((i + j) % 2 == 1) (A->elem)[i][j] *= (-1);
            }
        }
    }
    for (long long i = 0; i < n; ++i) {
        (f->elem)[0][i] = n1 * exp(M1 / (i + 1)) * cos(M1);
    }
}

void fill_matrix4 (matrix *A, matrix *f, long long n) {
    if (n != 100) {
        printf("ERROR ERROR ERROR\n");
        return;
    }
    long double M1 = 5, n1 = 100;
    long double qM = q_M(M1);
    for (long long i = 0; i < n; ++i) {
        for (long long j = 0; j < n; ++j) {
            if (i != j) {
                (A->elem)[i][j] = pow(qM, (double)(i + j + 2)) + 0.1 * ((double)j - (double)i);
            } else {
                (A->elem)[i][j] = pow(abs_d(qM - 1), (double)(i + j + 2));
                if ((i + j) % 2 == 1) (A->elem)[i][j] *= (-1);
            }
        }
    }
    for (long long i = 0; i < n; ++i) {
        (f->elem)[0][i] = abs_d(M1 - (n1 / 10)) * i * sin(M1);
    }
}

void fill_matrix5 (matrix *A, matrix *f, long long n) {
    if (n != 100) {
        printf("ERROR ERROR ERROR\n");
        return;
    }
    long double M1 = 6;
    long double qM = q_M(M1);
    for (long long i = 0; i < n; ++i) {
        for (long long j = 0; j < n; ++j) {
            if (i != j) {
                (A->elem)[i][j] = pow(qM, (double)(i + j + 2)) + 0.1 * ((double)j - (double)i);
            } else {
                (A->elem)[i][j] = pow(abs_d(qM - 1), (double)(i + j + 2));
                if ((i + j) % 2 == 1) (A->elem)[i][j] *= (-1);
            }
        }
    }
    for (long long i = 0; i < n; ++i) {
        (f->elem)[0][i] = M1 * exp(M1 / (i + 1)) * cos(M1 / (i + 1));
    }
}

void triangulate_matrix (matrix *matr, matrix *f) { //треугольный вид без выделения ведущего элемента
    long long m = matr->m;
    long long n = matr->n;
    long double **elem = matr->elem;
    for (long long i = 0; i < m; ++i) {
        for (long long j = i + 1; j < m; ++j) {
            long double koef = elem[i][j] / elem[i][i];
            for (long long k = i; k < n; ++k) {
                elem[k][j] -= elem[k][i] * koef;
            }
            if (f != NULL) {
                (f->elem)[0][j] -= (f->elem)[0][i] * koef;
            }
        }
    }
}

long double determinant (matrix *matr) { //поиск определителя
    long long n = matr->n;
    matrix *new = copy_matrix(matr);
    triangulate_matrix(new, NULL);
    long double res = 1.0;
    for (long long i = 0; i < n; ++i) {
        res *= (new->elem)[i][i];
    }
    delete_matrix(new);
    return res;
}

matrix *reverse_matrix (matrix *A) { //поиск обратной матрицы методом Гаусса
    long long m = A->m;
    long long n = A->n;
    if (m != n) return NULL;
    matrix *A1 = new_matrix(n, m);
    for (long long i = 0; i < n; ++i) {
        for (long long j = 0; j < m; ++j) {
            if (i == j) (A1->elem)[i][j] = 1; // инициализация
            else (A1->elem)[i][j] = 0;
        }
    }
    long double **elem = A->elem;
    long double **elem1 = A1->elem;
    for (long long i = 0; i < m; ++i) {
        for (long long j = 0; j < n; ++j) {
            elem1[j][i] /= elem[i][i];
        }
        for (long long j = 0; j < n; ++j) {
            if (j != i) elem[j][i] /= elem[i][i];
        }
        elem[i][i] = 1;
        for (long long j = i + 1; j < m; ++j) {
            long double koef = elem[i][j];
            for (long long k = 0; k < n; ++k) {
                elem[k][j] -= elem[k][i] * koef;
                elem1[k][j] -= elem1[k][i] * koef; // прямой ход
            }
        }
    }
    for (long long i = m - 1; i >= 0; --i) {
        for (long long j = i - 1; j >= 0; --j) {
            long double koef = elem[i][j];
            for (long long k = 0; k < n; ++k) {
                elem[k][j] -= elem[k][i] * koef;
                elem1[k][j] -= elem1[k][i] * koef; //обратный ход
            }
        }
    }
    return A1;
}

matrix *mult_matrix (matrix *A, matrix *B) { // перемножение матриц
    if (A->n != B->m) return NULL;
    matrix *res = new_matrix(A->m, B->n);
    long double **elem = res->elem;
    for (long long i = 0; i < B->n; ++i) { // n результата
        for (long long j = 0; j < A->m; ++j) { // m результата
            long double cur = 0;
            for (long long k = 0; k < A->n; ++k) {
                cur += (A->elem)[k][j] * (B->elem)[i][k];
            }
            elem[i][j] = cur;
        }
    }
    return res;
}

matrix *new_iter (matrix *x, matrix *A, matrix *f, long double omega) {
    matrix *A1 = copy_matrix(A);
    for (long long i = 0; i < A1->n; ++i) {
        for (long long j = 0; j <= i; ++j) {
            if (j == i) (A1->elem)[i][j] /= omega;
            else (A1->elem)[i][j] = 0;
        }
    }
    matrix *A1_rev = reverse_matrix(A1); // = (((1/omega) * D) + A(-)) ^ (-1)
    delete_matrix(A1);
    matrix *vect = mult_matrix(A, x);
    for (long long i = 0; i < vect->m; ++i) {
        (vect->elem)[0][i] = (f->elem)[0][i] - (vect->elem)[0][i]; // vect = (f - Ax_k)
    }
    matrix *res = mult_matrix(A1_rev, vect);
    delete_matrix(A1_rev);
    delete_matrix(vect);
    for (long long i = 0; i < res->m; ++i) {
        (res->elem)[0][i] += (x->elem)[0][i];
    }
    return res;
}

long double diff (matrix *f, matrix *x) { // среднеквадратичная ошибка
    long double err = 0;
    for (long long i = 0; i < f->m; ++i) {
        err += pow(((f->elem)[0][i] - (x->elem)[0][i]), 2);
    }
    return sqrt(err);
}
