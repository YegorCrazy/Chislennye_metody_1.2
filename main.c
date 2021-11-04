#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char **argv) {
    //long long n;
    //scanf("%lld", &n);
    unsigned n = 20;
    
    long double omega = 0.1;
    long double eps = 0.000001;
    
    matrix *A = new_matrix(n, n);
    matrix *f = new_matrix(n, 1);
    
    /*FILE *input = fopen("input.txt", "r");
    read_matrix(input, A);
    read_matrix(input, f);
    fclose(input);*/
    fill_matrix1(A, f, n);
    
    printf("Determinant is %Lf (if it is 0 or nan, result can't be calculated well)\n\n", determinant(A));
    
    while (omega < 2) {
        printf("Omega = %Lf\n", omega);
        matrix *x0 = new_matrix(n, 1);
        /*for (long long i = 0; i < x0->m; ++i) {
            (x0->elem)[0][i] = 1;
        }*/
        long long iter = 0;
        matrix *f1 = mult_matrix(A, x0);
        int flag = 0;
        
        while (abs_d(diff(f, f1)) > eps) {
            matrix *x1 = new_iter(x0, A, f, omega);
            delete_matrix(x0);
            x0 = x1;
            delete_matrix(f1);
            f1 = mult_matrix(A, x0);
            ++iter;
            if (iter > 100) {
                printf("Cycled\n\n");
                flag = 1;
                break;
            }
        }
        
        if (flag == 0) {
            printf("It took %Ld iterations to find the roots\n\n", iter);
            printf("The roots are:\n");
            print_matrix(x0);
            printf("The new right part is:\n");
            matrix *fin = mult_matrix(A, x0);
            print_matrix(fin);
            delete_matrix(fin);
        }
        
        delete_matrix(f1);
        delete_matrix(x0);
        
        omega += 0.1;
    }

    delete_matrix(A);
    delete_matrix(f);
    
    return 0;
}
