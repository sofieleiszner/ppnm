#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>


double Hilbert(int i, int j){
    double Hij = (1/(double)(i+j+1));
    return Hij;}

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++)printf("%10g ",gsl_vector_get(v,i));
	printf("\n");
}

void matrix_print(char s[], gsl_matrix* M){
    printf("%s\n",s);
	for(int i=0;i< M->size1;i++){
        for(int j = 0; j< M->size2; j++){
            //printf("i, j = %10d %10d \n", i,j);
        printf("%10g ", gsl_matrix_get(M,i,j));}
    printf("\n");}
    
}

int main(){
	int n=4;
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	//gsl_matrix* Acopy=gsl_matrix_alloc(n,n);
	for(int i=0; i< A->size1; i++)
		for(int j=0; j<A->size2; j++)
		{
		double Aij=Hilbert(i,j);
		gsl_matrix_set(A,i,j,Aij);
		}
	//gsl_matrix_memcpy(Acopy,A);
    matrix_print("The Hilbert matrix = ", A);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(n);
    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);
    gsl_eigen_symmv (A, eval, evec, w);
    //Sorterer dem, men det er nok ikke nødvendigt. 
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    //Printing the eigenvectors and values
    matrix_print("Eigenvectors", evec);
    vector_print("Eigenvalues", eval);   
    //Her kunne man evt. lave noget, som tjekker at Av=λv for alle egenvektorer

    
gsl_matrix_free(A);
//gsl_matrix_free(Acopy);
gsl_eigen_symmv_free(w);
return 0;
}