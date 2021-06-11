#include<gsl/gsl_matrix.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_vector.h>
#include <stdio.h>
#include<gsl/gsl_eigen.h>

void GenerateMatrix(gsl_matrix* A, int N) {
    
    for(int i = 0; i<N; i++) {
        for(int j = i; j<N; j++) {
			double tal=(double)rand()/RAND_MAX;
            gsl_matrix_set(A, i, j, tal);   
			gsl_matrix_set(A, j, i, tal);
		}
	}
    
}
void matrix_print(char s[], gsl_matrix* M){
    printf("%s\n",s);
	for(int i=0;i< M->size1;i++){
        for(int j = 0; j< M->size2; j++){
            //printf("i, j = %10d %10d \n", i,j);
        printf("%10g ", gsl_matrix_get(M,i,j));}
    printf("\n");}
    
}

void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int i=0;i<A->size1;i++){
		double new_aip=c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		double new_aiq=s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,new_aip);
		gsl_matrix_set(A,i,q,new_aiq);
		}
}


void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int j=0;j<A->size2;j++){
		double new_apj= c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
		double new_aqj=-s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,p,j,new_apj);
		gsl_matrix_set(A,q,j,new_aqj);
		}
}



void jacobi_diag(gsl_matrix* A, gsl_matrix* V){
int changed;
do{
	changed=0;
	for(int p=0;p<A->size1-1;p++)
	for(int q=p+1;q<A->size1;q++){
		double apq=gsl_matrix_get(A,p,q);
		double app=gsl_matrix_get(A,p,p);
		double aqq=gsl_matrix_get(A,q,q);
		double theta=0.5*atan2(2*apq,aqq-app);
		double c=cos(theta),s=sin(theta);
		double new_app=c*c*app-2*s*c*apq+s*s*aqq;
		double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
		if(new_app!=app || new_aqq!=aqq) // do rotation
			{
			changed=1;
			timesJ(A,p,q, theta);
			Jtimes(A,p,q,-theta); // A←J^T*A*J 
			timesJ(V,p,q, theta); // V←V*J
			}
	}
}while(changed!=0);
}


void JtAJ(gsl_matrix* A, int p, int q, double theta) {
   double c = cos(theta);
   double s = sin(theta);

   double apq = gsl_matrix_get(A, p, q);
   double app = gsl_matrix_get(A, p, p);
   double aqq = gsl_matrix_get(A, q, q);
   gsl_matrix_set(A, p, q, 0);
   gsl_matrix_set(A, p, p, c*c*app - 2*s*c*apq + s*s*aqq);
   gsl_matrix_set(A, q, q, s*s*app + 2*s*c*apq + c*c*aqq);

   // loop over remaining upper triangle elements in columns p, q and rows p, q
   for (int i = 0; i < p; ++i) {
      double aip = gsl_matrix_get(A, i, p);
      double aiq = gsl_matrix_get(A, i, q);
      gsl_matrix_set(A, i, p, c*aip - s*aiq);
      gsl_matrix_set(A, i, q, s*aip + c*aiq);
   }

   for (int i = p + 1; i < q; ++i) {
      double api = gsl_matrix_get(A, p, i);
      double aiq = gsl_matrix_get(A, i, q);
      gsl_matrix_set(A, p, i, c*api - s*aiq);
      gsl_matrix_set(A, i, q, s*api + c*aiq);
   }

   for (int i = q + 1; i < A->size1; ++i) {
      double api = gsl_matrix_get(A, p, i);
      double aqi = gsl_matrix_get(A, q, i);
      gsl_matrix_set(A, p, i, c*api - s*aqi);
      gsl_matrix_set(A, q, i, s*api + c*aqi);
   }

}


void jacobi_diag_op(gsl_matrix* A, gsl_matrix* V){
int changed;
do{
	changed=0;
	for(int p=0;p<A->size1-1;p++)
	for(int q=p+1;q<A->size1;q++){
		double apq=gsl_matrix_get(A,p,q);
		double app=gsl_matrix_get(A,p,p);
		double aqq=gsl_matrix_get(A,q,q);
		double theta=0.5*atan2(2*apq,aqq-app);
		double c=cos(theta),s=sin(theta);
		double new_app=c*c*app-2*s*c*apq+s*s*aqq;
		double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
		if(new_app!=app || new_aqq!=aqq) // do rotation
			{
			changed=1;
			JtAJ(A, p, q, theta); // A <- Jt*A*J
			timesJ(V,p,q, theta); // V←V*J
			}
		
	}

}while(changed!=0);
}



void builtH(gsl_matrix* H){
	int n=H->size1;
	double s=1.0/(n+1);

	for(int i=0;i<n-1;i++){
		gsl_matrix_set(H,i,i,-2);
	
	gsl_matrix_set(H,i,i+1,1);
	gsl_matrix_set(H,i+1,i,1);
	}
	gsl_matrix_set(H,n-1,n-1,-2);
	gsl_matrix_scale(H,-1/s/s);

}

double exact(double x, int n){
	double value=sqrt(2)*0.145*sin(n*M_PI*x);
	return value;
}


void MeasureTime(int Nmin,int Nmax) {
    FILE * f3=fopen("out_time_jacobi.txt","w");    
    FILE * f4=fopen("out_time_gsl.txt","w");
	FILE * f5=fopen("out_time_jacobiop.txt","w");
	
    for (int N=Nmin; N<Nmax; N += 5){

		gsl_matrix* A = gsl_matrix_alloc(N, N); 
		gsl_matrix* A1 = gsl_matrix_alloc(N+1, N+1); 
		gsl_matrix* A2 = gsl_matrix_alloc(N+2, N+2); 
		gsl_matrix* A3 = gsl_matrix_alloc(N+3, N+3); 
		gsl_matrix* A4 = gsl_matrix_alloc(N+4, N+4); 
		
		gsl_matrix* V = gsl_matrix_alloc(N, N); 
		gsl_matrix* V1 = gsl_matrix_alloc(N+1, N+1); 
		gsl_matrix* V2 = gsl_matrix_alloc(N+2, N+2); 
		gsl_matrix* V3 = gsl_matrix_alloc(N+3, N+3); 
		gsl_matrix* V4 = gsl_matrix_alloc(N+4, N+4); 
      
        GenerateMatrix(A, N);
        GenerateMatrix(A1, N+1);
        GenerateMatrix(A2, N+2);
        GenerateMatrix(A3, N+3);
        GenerateMatrix(A4, N+4);
			
     	gsl_matrix_set_identity(V); 
     	gsl_matrix_set_identity(V1); 
     	gsl_matrix_set_identity(V2); 
     	gsl_matrix_set_identity(V3); 
     	gsl_matrix_set_identity(V4); 
        
		clock_t begin = clock(); 
        jacobi_diag(A,V);        
        jacobi_diag(A1,V1);        
        jacobi_diag(A2,V2);        
        jacobi_diag(A3,V3);        
        jacobi_diag(A4,V4);        
		clock_t end = clock();
        fprintf(f3, "%3d  %Lf \n", N, (long double)(end - begin) / CLOCKS_PER_SEC);

		
		gsl_matrix* W = gsl_matrix_alloc(N, N); 
		gsl_matrix* W1 = gsl_matrix_alloc(N+1, N+1); 
		gsl_matrix* W2 = gsl_matrix_alloc(N+2, N+2); 
		gsl_matrix* W3 = gsl_matrix_alloc(N+3, N+3); 
		gsl_matrix* W4 = gsl_matrix_alloc(N+4, N+4); 

		gsl_vector *S=gsl_vector_alloc(N);
		gsl_vector *S1=gsl_vector_alloc(N+1);
		gsl_vector *S2=gsl_vector_alloc(N+2);
		gsl_vector *S3=gsl_vector_alloc(N+3);
		gsl_vector *S4=gsl_vector_alloc(N+4);
		
		gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(N);
		gsl_eigen_symmv_workspace* w1 = gsl_eigen_symmv_alloc(N+1);
		gsl_eigen_symmv_workspace* w2 = gsl_eigen_symmv_alloc(N+2);
		gsl_eigen_symmv_workspace* w3 = gsl_eigen_symmv_alloc(N+3);
		gsl_eigen_symmv_workspace* w4 = gsl_eigen_symmv_alloc(N+4);

        
		GenerateMatrix(A, N);
        GenerateMatrix(A1, N+1);
        GenerateMatrix(A2, N+2);
        GenerateMatrix(A3, N+3);
        GenerateMatrix(A4, N+4); 


        clock_t begin_gsl = clock(); 
        gsl_eigen_symmv(A, S, W, w);
        gsl_eigen_symmv(A1, S1, W1, w1);
        gsl_eigen_symmv(A2, S2, W2, w2);
        gsl_eigen_symmv(A3, S3, W3, w3);
        gsl_eigen_symmv(A4, S4, W4, w4);
		clock_t end_gsl = clock();
        fprintf(f4, "%3d  %Lf \n", N, (long double)(end_gsl - begin_gsl) / CLOCKS_PER_SEC);
        

		GenerateMatrix(A, N);
        GenerateMatrix(A1, N+1);
        GenerateMatrix(A2, N+2);
        GenerateMatrix(A3, N+3);
        GenerateMatrix(A4, N+4);
			
     	gsl_matrix_set_identity(V); 
     	gsl_matrix_set_identity(V1); 
     	gsl_matrix_set_identity(V2); 
     	gsl_matrix_set_identity(V3); 
     	gsl_matrix_set_identity(V4); 
		
		clock_t begin_op = clock(); 
        jacobi_diag_op(A,V);        
        jacobi_diag_op(A1,V1);        
        jacobi_diag_op(A2,V2);        
        jacobi_diag_op(A3,V3);        
        jacobi_diag_op(A4,V4);        
		clock_t end_op = clock();
        
        fprintf(f5, "%3d  %Lf \n", N, (long double)(end_op - begin_op) / CLOCKS_PER_SEC);

		
        gsl_matrix_free(A);
        gsl_matrix_free(A1);
        gsl_matrix_free(A2);
        gsl_matrix_free(A3);
        gsl_matrix_free(A4);
		
        gsl_matrix_free(V);
        gsl_matrix_free(V1);
        gsl_matrix_free(V2);
        gsl_matrix_free(V3);
        gsl_matrix_free(V4);
		
        gsl_matrix_free(W);
        gsl_matrix_free(W1);
        gsl_matrix_free(W2);
        gsl_matrix_free(W3);
        gsl_matrix_free(W4);
		
        gsl_vector_free(S);
        gsl_vector_free(S1);
        gsl_vector_free(S2);
        gsl_vector_free(S3);
        gsl_vector_free(S4);
		
    
		gsl_eigen_symmv_free(w);
		gsl_eigen_symmv_free(w1);
		gsl_eigen_symmv_free(w2);
		gsl_eigen_symmv_free(w3);
		gsl_eigen_symmv_free(w4);
		}
	
	fclose(f3);
	fclose(f4);


}