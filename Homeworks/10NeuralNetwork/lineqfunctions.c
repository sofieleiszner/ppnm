#include <math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_linalg.h>


//General functions: 

void matrix_print(char s[], gsl_matrix* M){
    printf("%s\n",s);
	for(int i=0;i< M->size1;i++){
        for(int j = 0; j< M->size2; j++){
            //printf("i, j = %10d %10d \n", i,j);
        printf("%10g ", gsl_matrix_get(M,i,j));}
    printf("\n");}
    
}


void GenerateMatrix(gsl_matrix* A, int N, int M) {
    
    for(int i = 0; i<N; i++) {
        for(int j = 0; j<M; j++) {
            gsl_matrix_set(A, i, j, (double)rand()/RAND_MAX);}   
    }
    
}

//GS_decomp


void GS_decomp(gsl_matrix* A, gsl_matrix* R) {
   int N = A->size1; //Size1 er antal rækker
   int M = A->size2;   //Size1 er antal søjler
   for (int i =0; i<M; i++){
       
       //Sets Rii=sqrt(ai^T*ai)
       gsl_vector* ai = gsl_vector_alloc(N); 
       gsl_matrix_get_col(ai, A, i);
       double dot; 
       gsl_blas_ddot(ai, ai, &dot); 
       gsl_matrix_set(R, i, i, sqrt(dot));
       
       //Sets qi = ai/Rii
       gsl_vector_scale(ai, (double)1/sqrt(dot)); 
       gsl_matrix_set_col(A, i, ai);
       
       
       for (int j = i+1; j < M; j++){ 
            
            //Sets Rij = qi^T*aj
            double dot2; 
            gsl_vector* aj = gsl_vector_alloc(N); 
            gsl_matrix_get_col(aj, A, j);
            //gsl_matrix_get_col(ai, A, i);
            gsl_blas_ddot(ai, aj, &dot2); 
            gsl_matrix_set(R, i, j, dot2);
            
            //Sets aj = aj-qi*Rij
            gsl_vector_scale(ai,dot2);
            gsl_vector_sub(aj, ai); 
            gsl_matrix_set_col(A, j, aj); //aj = aj-qi*Rij
            gsl_vector_scale(ai,(double)1/dot2);
            gsl_vector_free(aj);
            
            }
    gsl_vector_free(ai);    
   }


}



//GS_solve
// i is the rows, j is the coloumns.

double Rsum(gsl_matrix *R, gsl_vector *x, int i){
	if(i==x->size-1) return 0;
    
	double sum=0;
	int j=x->size-1; //-1 to ensure proper indexing
	for(j=j;j>i;j--)
		{
		sum+=gsl_matrix_get(R,i,j)*gsl_vector_get(x,j);
		}
return sum;
}


void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){

	gsl_matrix_transpose(Q);
	gsl_blas_dgemv(CblasNoTrans, 1.0, Q, b, 0.0, x);

	int i=x->size-1; //size-1 to ensure proper indexing. i is the rows in R.
	for(i=i;i>=0;i--)
		{
		gsl_vector_set(x,i,(gsl_vector_get(x,i)-Rsum(R,x,i))/gsl_matrix_get(R,i,i));
		}
    gsl_matrix_transpose(Q);

//	xn=yn/Rnn
//	xn1=(yn1- Rn1n xn)/Rn1n1
//	xn2=(yn2- Rn2n xn -Rn2n1 xn1)/Rn2n2
//	xn3=(yn3- Rn3n xn -Rn3n1 xn1 -Rn3n2 xn2)/Rn3n3
//	xn4=(yn4- Rn4n xn -Rn4n1 xn1 -Rn4n2 xn2 -Rn4n3 xn3)/Rn4n4
}


//GS_inverse  
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
	gsl_matrix_set_identity(B); 
	for(int i=0; i<Q->size2; i++){
		gsl_vector_view v = gsl_matrix_column(B,i);
        gsl_vector* v_save = gsl_vector_alloc(B->size1);
        gsl_vector_memcpy(v_save,&v.vector);
		GS_solve(Q, R, v_save, &v.vector);
        gsl_vector_free(v_save);
	}	
    
}



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
        
        fprintf(f1, "%3d %Lf \n", N, (long double)(end - begin) / CLOCKS_PER_SEC);

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
        fprintf(f2, "%3d %Lf \n", N, (long double)(end_gsl - begin_gsl) / CLOCKS_PER_SEC);

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


