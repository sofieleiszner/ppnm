#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <gsl/gsl_blas.h>


void GenerateMatrix(gsl_matrix* A, int N, int M);
void matrix_print(char s[], gsl_matrix* M);
void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void MeasureDecompTime(int Nmax);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void vector_print(char s[], gsl_vector* v); 
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);

int main(){
    
    int N = 6; int M = 6; 

    gsl_matrix* A = gsl_matrix_alloc(N, M); 
    gsl_matrix* A_save = gsl_matrix_alloc(N, M); 
    gsl_matrix* R = gsl_matrix_alloc(M, M); 
    

    GenerateMatrix(A, N, M);
    gsl_matrix_memcpy(A_save,A);
    printf("TASK A.1 \n"); 
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

    gsl_vector* b=gsl_vector_alloc(N);
	for (int i=0;i<N;i++)
		gsl_vector_set(b,i,(double)rand()/RAND_MAX);
    gsl_vector* x = gsl_vector_alloc(N);
    
    printf("\nTASK A.2 \n"); 

    GS_solve(A, R, b, x);   
    //Checks Ax = b   
    gsl_vector* Ax = gsl_vector_alloc(N);
    gsl_blas_dgemv(CblasNoTrans, 1.0, A_save, x, 0.0, Ax); 
    matrix_print("A = ", A_save); 
    vector_print("x=", x);
    vector_print("b=", b);
    printf("Checks that Ax = b:\n"); 
    vector_print("Ax =", Ax); 
    
    printf("\nTASK B \n");
    gsl_matrix* B = gsl_matrix_alloc(N,N);
    GS_inverse(A, R, B);
    matrix_print("Inverse of A, B=", B);

    printf("Checks that AB = I"); 
    gsl_matrix* AB = gsl_matrix_alloc(N,N); 
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A_save, B, 0.0, AB); 
    matrix_print("AB =", AB); 
    printf("Checks that BA = I"); 
    gsl_matrix* BA = gsl_matrix_alloc(N,N); 
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, B, A_save, 0.0, BA); 
    matrix_print("BA =", BA); 

    printf("\nTASK C \n");
    printf("For time measuring, see .png file \n"); 
    
    
    gsl_matrix_free(A); 
    gsl_matrix_free(R);
    gsl_matrix_free(QTQ);
    gsl_matrix_free(AT);
    gsl_matrix_free(QR);
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_vector_free(Ax);
    gsl_matrix_free(B);
    gsl_matrix_free(AB);
    gsl_matrix_free(A_save);
    
    int Nmax = 500; 
    MeasureDecompTime(Nmax);
    
return 0; 
}