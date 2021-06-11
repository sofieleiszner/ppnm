#include <math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

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