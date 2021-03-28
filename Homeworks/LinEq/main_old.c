#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <time.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>


void GenerateMatrix(gsl_matrix* A, int N, int M);
void matrix_print(char s[], gsl_matrix* M);
void GS_decomp(gsl_matrix* A, gsl_matrix* R);

int main() {
    int Nmax = 200;
    clock_t begin = clock();
    for (int N = 2; N<Nmax; N++){
    gsl_matrix* A = gsl_matrix_alloc(N, N); 
    gsl_matrix* R = gsl_matrix_alloc(N, N); 
    GenerateMatrix(A, N, N);    
    GS_decomp(A, R);
    clock_t end = clock();
    
    printf("%3d  %Lf \n", N, (long double)(end - begin) / CLOCKS_PER_SEC);
    }
    FILE * f=fopen("out_gsl.txt","w");
    clock_t begin_gsl = clock();
    for (int N = 2; N<Nmax; N++){
    gsl_matrix* A = gsl_matrix_alloc(N, N); 
    GenerateMatrix(A, N, N);
    gsl_vector* tau = gsl_vector_alloc(N);
    gsl_linalg_QR_decomp(A, tau);
    clock_t end_gsl = clock();
    
    fprintf(f, "%3d  %Lf \n", N, (long double)(end_gsl - begin_gsl) / CLOCKS_PER_SEC);
    }
    
    //gsl_matrix_free(A); 
    //gsl_matrix_free(R);

    
return 0; 
}