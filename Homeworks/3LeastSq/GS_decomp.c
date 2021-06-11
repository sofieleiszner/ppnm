#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_blas.h>

void matrix_print(char s[], gsl_matrix* M){
    printf("%s\n",s);
	for(int i=0;i< M->size1;i++){
        for(int j = 0; j< M->size2; j++){
            //printf("i, j = %10d %10d \n", i,j);
        printf("%10g ", gsl_matrix_get(M,i,j));}
    printf("\n");}
    
}

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++)printf("%10g \n",gsl_vector_get(v,i));
//	printf("\n");
}

void GenerateMatrix(gsl_matrix* A, int N, int M) {
    
    for(int i = 0; i<N; i++) {
        for(int j = 0; j<M; j++) {
            gsl_matrix_set(A, i, j, (double)rand()/RAND_MAX);}   
    }
    
}



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