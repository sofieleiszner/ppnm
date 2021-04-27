#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <gsl/gsl_linalg.h>

#include <math.h>
#include <gsl/gsl_blas.h>

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void GS_decomp(gsl_matrix* A, gsl_matrix* R);



void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps){
    int n = x->size;
    gsl_matrix* J = gsl_matrix_alloc(n,n);
    gsl_vector* fx = gsl_vector_alloc(n);
    gsl_vector* dfx = gsl_vector_alloc(n);
    gsl_vector* df = gsl_vector_alloc(n);     
    gsl_vector* fy = gsl_vector_alloc(n);     
    gsl_vector* Dx = gsl_vector_alloc(n);    
    gsl_vector* y = gsl_vector_alloc(n);    
    gsl_matrix* R = gsl_matrix_alloc(n,n);
    //gsl_matrix* QR = gsl_matrix_alloc(n,n);
    double lambda; 
    double dx = 1e-6;
    
    
    while(1){

        f(x, fx); 
        for (int j =0; j<n; j++) {
            gsl_vector_set(x,j,gsl_vector_get(x,j)+dx);
            f(x,dfx);
            gsl_vector_memcpy(df,dfx);
            gsl_vector_sub(df,fx);
            for (int i=0; i<n; i++){
                gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/dx);             
            }
        gsl_vector_set(x,j,gsl_vector_get(x,j)-dx);   
        }
        
        
        GS_decomp(J,R);
        gsl_vector_memcpy(Dx,x);
        gsl_vector_scale(fx,-1);
        GS_solve(J,R,fx,Dx);
        gsl_vector_scale(fx,-1);
    
        lambda = 2;
        while(1){
            lambda = (double)lambda/2;
            gsl_vector_scale(Dx,lambda);
            gsl_vector_memcpy(y,x); 
            gsl_vector_add(y,Dx);
            f(y,fy);
            gsl_vector_scale(Dx,(double)1/lambda);
        
            if(gsl_blas_dnrm2(fy)<(1-(double)lambda/2)*gsl_blas_dnrm2(fx) || lambda < 0.02) {
                break;
            }
        }
        gsl_vector_memcpy(x,y); 
        gsl_vector_memcpy(fx,fy); 
        if(gsl_blas_dnrm2(Dx)<dx || gsl_blas_dnrm2(fx) < eps) {
                break;
                }
    }
    
    
    gsl_matrix_free(J  );
    gsl_vector_free(fx );
    gsl_vector_free(dfx);
    gsl_vector_free(df ); 
    gsl_vector_free(fy ); 
    gsl_vector_free(Dx );
    gsl_vector_free(y  );
    gsl_matrix_free(R);
    
    
}


