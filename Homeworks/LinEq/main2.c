#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <gsl/gsl_blas.h>


void GenerateMatrix(gsl_matrix* A, int N, int M);
void matrix_print(char s[], gsl_matrix* M);
void GS_decomp(gsl_matrix* A, gsl_matrix* R);

int main() {
    int N = 5; int M = 3; 

    gsl_matrix* A = gsl_matrix_alloc(N, M); 
    gsl_matrix* R = gsl_matrix_alloc(M, M); 
    
    GenerateMatrix(A, N, M);
    matrix_print("Random matrix A = ", A);
    
    GS_decomp(A, R);

    printf("\nResults after decomposition: \n");
    matrix_print("R = ", R);
    matrix_print("Q = ", A);
    
    printf("\nChecks:\n");
    //check that Q^T*Q;
    printf("Check that Q^T*Q = 1: \n");
    // gsl_blas_dgemm()
    gsl_matrix* QTQ = gsl_matrix_alloc(M, M); 
    gsl_matrix* AT = gsl_matrix_alloc(M, N);
    gsl_matrix_transpose_memcpy(AT, A); 
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AT, A, 0.0, QTQ);
    matrix_print("Q^T*Q =", QTQ); 
    

    printf("Check that QR=A:\n");
    gsl_matrix* QR = gsl_matrix_alloc(N, M); 
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, R, 0.0, QR); 
    matrix_print("QR=", QR);
    
    gsl_matrix_free(A); 
    gsl_matrix_free(R);
    gsl_matrix_free(QTQ);
    gsl_matrix_free(AT);
    gsl_matrix_free(QR);
    
return 0; 
}