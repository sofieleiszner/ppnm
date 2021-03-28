#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>


void GenerateMatrix(gsl_matrix* A, int N, int M);
void GS_decomp(gsl_matrix* A, gsl_matrix* R);


void MeasureDecompTime(int Nmax) {
    FILE * f1=fopen("out_decomp.txt","w");    
    FILE * f2=fopen("out_gsl.txt","w");

    for (int N = 0; N<Nmax; N += 5){
        
        gsl_matrix* A = gsl_matrix_alloc(N, N); 
        gsl_matrix* A1 = gsl_matrix_alloc(N+1, N+1); 
        gsl_matrix* A2 = gsl_matrix_alloc(N+2, N+2); 
        gsl_matrix* A3 = gsl_matrix_alloc(N+3, N+3); 
        gsl_matrix* A4 = gsl_matrix_alloc(N+4, N+4); 

        gsl_matrix* R = gsl_matrix_alloc(N, N); 
        gsl_matrix* R1 = gsl_matrix_alloc(N+1, N+1); 
        gsl_matrix* R2 = gsl_matrix_alloc(N+2, N+2); 
        gsl_matrix* R3 = gsl_matrix_alloc(N+3, N+3); 
        gsl_matrix* R4 = gsl_matrix_alloc(N+4, N+4); 
        
        GenerateMatrix(A , N, N);   
        GenerateMatrix(A1, N+1, N+1);   
        GenerateMatrix(A2, N+2, N+2);   
        GenerateMatrix(A3, N+3, N+3);   
        GenerateMatrix(A4, N+4, N+4); 
        
        clock_t begin = clock(); 
        GS_decomp(A, R);        
        GS_decomp(A1, R1);   
        GS_decomp(A2, R2);   
        GS_decomp(A3, R3);   
        GS_decomp(A4, R4);   
        clock_t end = clock();
        
        fprintf(f1, "%3d  %Lf \n", N, (long double)(end - begin) / CLOCKS_PER_SEC);

        gsl_matrix_free(R); 
        gsl_matrix_free(R1); 
        gsl_matrix_free(R2); 
        gsl_matrix_free(R3); 
        gsl_matrix_free(R4); 

        gsl_vector* tau = gsl_vector_alloc(N);
        gsl_vector* tau1 = gsl_vector_alloc(N+1);
        gsl_vector* tau2 = gsl_vector_alloc(N+2);
        gsl_vector* tau3 = gsl_vector_alloc(N+3);
        gsl_vector* tau4 = gsl_vector_alloc(N+4);

        clock_t begin_gsl = clock(); 
        
        gsl_linalg_QR_decomp(A, tau);
        gsl_linalg_QR_decomp(A1, tau1);
        gsl_linalg_QR_decomp(A2, tau2);
        gsl_linalg_QR_decomp(A3, tau3);
        gsl_linalg_QR_decomp(A4, tau4);

        clock_t end_gsl = clock();
        fprintf(f2, "%3d  %Lf \n", N, (long double)(end_gsl - begin_gsl) / CLOCKS_PER_SEC);

        gsl_matrix_free(A);
        gsl_matrix_free(A1);
        gsl_matrix_free(A2);
        gsl_matrix_free(A3);
        gsl_matrix_free(A4);
                
        gsl_vector_free(tau);
        gsl_vector_free(tau1);
        gsl_vector_free(tau2);
        gsl_vector_free(tau3);
        gsl_vector_free(tau4);
        }
        


}