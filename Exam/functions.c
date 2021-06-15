#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<stdarg.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>


void GenerateMatrix(gsl_matrix* A, int M, int N) {
    
    for(int i = 0; i<M; i++) {
        for(int j = 0; j<N; j++) {
            gsl_matrix_set(A, i, j, (double)rand()/RAND_MAX);}   
    }
    
}

void GenerateVector(gsl_vector* v){
	for (int i=0;i<v->size;i++)
		gsl_vector_set(v,i,(double)rand()/RAND_MAX);
}

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++)printf("%10g \n",gsl_vector_get(v,i));
//	printf("\n");
}

void matrix_print(char s[], gsl_matrix* M){
    printf("%s\n",s);
	for(int i=0;i< M->size1;i++){
        for(int j = 0; j< M->size2; j++){
            if (fabs(gsl_matrix_get(M,i,j))<pow(10,-8)){
                printf("%6f ", 0.000000000);}
            else {
                printf("%6f ", gsl_matrix_get(M,i,j));}
        }
        printf("\n");}

        }


double vnorm(gsl_vector* v){ //Takes gsl_vector and calculates its norm (length)
    double sum = 0;
    for (int i=0;i<v->size;i++){
        sum += gsl_vector_get(v, i)*gsl_vector_get(v, i); 
    }
    return sqrt(sum); 
}



void GKLBidiag(gsl_matrix* A, gsl_matrix* U, gsl_matrix* V, gsl_matrix* B){
 
    int M = A->size1;   //M er antal rækker
    int N = A->size2; //N er antal søjler
    
   
    //Step 1: Sets v1 = unit 2-norm vector
    
    gsl_vector* vk = gsl_vector_alloc(N); 
    gsl_vector_set(vk,0,1); 
    for(int i = 1; i<vk->size; i++){ 
        gsl_vector_set(vk,i,0); 
        }
    gsl_matrix_set_col(V, 0, vk); 
    
    //Step 1: Sets β0 = 0. Needs to set these into B in the end a long with α
    gsl_vector* beta = gsl_vector_alloc(N+1); //A vector of the beta-values
    gsl_vector* alpha = gsl_vector_alloc(N); //A vector of the beta-values
    gsl_vector_set(beta,0,0); 
    
    
    //Allocates before loop
    gsl_vector* Avk = gsl_vector_alloc(M);
    gsl_vector* Atuk = gsl_vector_alloc(N);
    double dot; 
    gsl_vector* ai = gsl_vector_alloc(N);
    gsl_vector* aj = gsl_vector_alloc(M);
    gsl_vector* uk = gsl_vector_alloc(M);
    gsl_matrix* At = gsl_matrix_alloc(N,M); 
    gsl_matrix_transpose_memcpy(At,A); 

    //Sets uk to zero and puts it into U's first column
    gsl_vector_set_zero(uk); 
    gsl_matrix_set_col(U,0,uk); 


    // Step 2-8: 
    for (int k = 1; k<N+1; k++) { //    

        //Step 3: Sets uk = A*v_k-beta_k-1*u_k-1
        gsl_matrix_get_col(vk, V, k-1);  //vk

        //Calculates A*vk: 
        for (int i =0; i<M; i++){
            gsl_matrix_get_row(ai,A,i); // Takes i'th row of A into ai
            gsl_blas_ddot(ai, vk, &dot); //Puts dot product of ai and ck into dot
            gsl_vector_set(Avk, i, dot);  //Puts it into vector Ack
        }

        //gsl_vector_axpby(1, Avk, gsl_vector_get(beta, k-1), uk); //Puts result into 4th argument. Arguments: α, x, β, y. y <- αx+βy = 1*Avk+β_k-1*u_k-1
        if (k==1) {gsl_matrix_get_col(uk, U, k-1);}
        else {gsl_matrix_get_col(uk, U, k-2);}
        gsl_blas_daxpy(-gsl_vector_get(beta, k-1), uk, Avk); //
        gsl_vector_memcpy(uk, Avk); 

        //Step 4: αk = ||uk|| 
        gsl_vector_set(alpha, k-1, vnorm(uk));
        
        //Step 5: uk = uk/αk
        gsl_vector_scale(uk, 1/gsl_vector_get(alpha, k-1)); 
        
        //Step 6: vk+1 = At*uk-αk*vk
        
        //At*uk is calulated and put into Auk
        for (int i =0; i<N; i++){
            gsl_matrix_get_row(aj,At,i); 
            gsl_blas_ddot(aj, uk, &dot); 
            gsl_vector_set(Atuk, i, dot);
        }
        
        //Wants to calculate -α_k*v_k+Atuk
        gsl_blas_daxpy(-gsl_vector_get(alpha, k-1), vk, Atuk); //Puts result into 4th argument. Arguments: α, x, β, y. y <- αx+βy = 1*Avk+β_k-1*u_k-1
        gsl_vector_memcpy(vk, Atuk);
        

        //Step 7: βk = ||vk+1||
        gsl_vector_set(beta,k,vnorm(vk)); 
        
        //Step 8: vk+1 = vk+1/βk
        gsl_vector_scale(vk, 1/gsl_vector_get(beta,k)); 
        
        gsl_matrix_set(B,k-1,k-1,gsl_vector_get(alpha,k-1)); 
        
        gsl_matrix_set_col(U, k-1, uk);  
        if (k<N){
            gsl_matrix_set_col(V, k, vk);  
            gsl_matrix_set(B,k-1,k,gsl_vector_get(beta,k));
            }
    }
    
    //Frees vectors and matrices
    gsl_vector_free(beta);  
    gsl_vector_free(alpha); 
    gsl_vector_free(Avk );
    gsl_vector_free(vk  );
    gsl_vector_free(Atuk);
    gsl_vector_free(ai  );
    gsl_vector_free(aj  );
    gsl_vector_free(uk  );
    gsl_matrix_free(At  ); 
}


double detfromB(gsl_matrix* B){ //B is a upper triangular matrix. 
    //The determinant of A can be found as the product of diagonal elements of B.
    assert(B->size1==B->size2); 
    double sum = 1.0; 
    for (int i =0; i<B->size1; i++){
        sum*=gsl_matrix_get(B,i,i);     
    }
    return sum; 
}



double detfromA(gsl_matrix* A){  
    //The determinant of A can be found as the product of diagonal elements of B. 
    
    assert(A->size1==A->size2); 

    int M = A->size1;   //M er antal rækker
    int N = A->size2; //N er antal søjler
    
    gsl_matrix* V = gsl_matrix_alloc(N,N); 
    gsl_matrix* U = gsl_matrix_alloc(M,M); 
    gsl_matrix* B = gsl_matrix_alloc(M,N); //M,N 
        
    GKLBidiag(A, U, V, B); 
     
    double sum = 1.0; 
    for (int i =0; i<B->size1; i++){
        sum*=gsl_matrix_get(B,i,i);     
    }
    
    gsl_matrix_free(V); 
    gsl_matrix_free(U); 
    gsl_matrix_free(B); 
    return sum; 
}


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


void SolveGS(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){ //Solves QRx=b

	gsl_blas_dgemv(CblasTrans, 1.0, Q, b, 0.0, x);

	int i=x->size-1; //size-1 to ensure proper indexing. i is the rows in R.
	for(i=i;i>=0;i--)
		{
		gsl_vector_set(x,i,(gsl_vector_get(x,i)-Rsum(R,x,i))/gsl_matrix_get(R,i,i));
		}

}

void SolveGKL(gsl_matrix* U, gsl_matrix* B, gsl_matrix* V, gsl_vector* b, gsl_vector* y){
    //Wants to solve UBV^Tx=b
    //B*V^T*x=U^T*b
    //Trick: Solves B*y = c for y, where V^T*x=y and c = U^T*b
    SolveGS(U, B, b, y); 
    ////Since V^T*x = y, x = V*y
    gsl_vector* Vy = gsl_vector_alloc(y->size); 
    gsl_blas_dgemv(CblasNoTrans, 1.0, V, y, 0.0, Vy); 
    gsl_vector_memcpy(y,Vy); 
    gsl_vector_free(Vy); 
}

void InversGKL(gsl_matrix* U, gsl_matrix* B, gsl_matrix* V, gsl_matrix* Ainv){
	gsl_matrix_set_identity(Ainv); 
	for(int i=0; i<U->size2; i++){
		gsl_vector_view v = gsl_matrix_column(Ainv,i);
        gsl_vector* v_save = gsl_vector_alloc(Ainv->size1);
        gsl_vector_memcpy(v_save,&v.vector);
		SolveGKL(U, B, V, v_save, &v.vector);
        gsl_vector_free(v_save);
	}	
    
}

