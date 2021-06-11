#include <math.h>
#include <stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

void matrix_print(char s[], gsl_matrix* M);
void vector_print(char s[], gsl_vector* v); 
void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);

double f(int k, double z){
	if (k == 0) {return 1.0;}
	if (k == 1) {return -z;}

return 0;
}

void leastsq(gsl_vector* x, gsl_vector* y, gsl_vector* dy, gsl_vector* c,gsl_vector* dc, double f(int k, double z)){
	
	gsl_matrix* A = gsl_matrix_alloc(x->size,c->size);
	gsl_matrix* R = gsl_matrix_alloc(c->size,c->size);
	
	for (int i=0; i<A->size1; i++){
		for (int k=0; k<A->size2; k++){
			gsl_matrix_set(A, i, k, f(k, gsl_vector_get(x,i))/gsl_vector_get(dy,i));
		}
	}
	
	matrix_print("A=",A); 
	
	gsl_vector* b = gsl_vector_alloc(dy->size);
	gsl_vector_memcpy(b,y);
	gsl_vector_div(b, dy);
	
	GS_decomp(A, R);
	GS_solve(A, R, b, c);
	
	gsl_matrix* Sigma = gsl_matrix_alloc(c->size,c->size);
	gsl_matrix* Rsigma = gsl_matrix_alloc(c->size,c->size);
	gsl_matrix* Sigma_inverse = gsl_matrix_alloc(c->size,c->size);
	
	//matrix_print("R = ", R);
	
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, R, R, 0.0, Sigma);
	//matrix_print("Sigma =", Sigma);
	
	GS_decomp(Sigma, Rsigma);
	//matrix_print("Sigma_decomp =", Sigma);
	GS_inverse(Sigma,Rsigma,Sigma_inverse);
	matrix_print("Calculated Covariance marix =", Sigma_inverse);
	for (int i=0;i<c->size;i++){ //bestemmer usikkerhederne på a og b
		gsl_vector_set(dc,i,sqrt(gsl_matrix_get(Sigma_inverse,i,i)));
	}
	
	
	
	
}