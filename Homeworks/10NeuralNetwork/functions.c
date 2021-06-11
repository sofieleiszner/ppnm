#include <math.h>
#include <float.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>

static const double DELTA=1.0/524288;
void numeric_gradient
(double F(gsl_vector*), gsl_vector*x, gsl_vector*grad){
	double fx=F(x);
	for(int i=0;i<x->size;i++){
		double dx,xi=gsl_vector_get(x,i);
		if(fabs(xi)<sqrt(DELTA)) dx=DELTA;
		else dx=fabs(xi)*DELTA;
		gsl_vector_set(x,i,xi+dx);
		gsl_vector_set(grad,i,(F(x)-fx)/dx);
		gsl_vector_set(x,i,xi);
	}
}

void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void matrix_print(char s[], gsl_matrix* M);
void vector_print(char s[], gsl_vector* v);


void qnewton(
	double f(gsl_vector* x), // objective function 
	gsl_vector* x, //  on input: starting point, on exit: approximation to root 
	double eps   //accuracy goal, on exit |gradient| should be <eps 
	){
	//vector_print("x i qnewton",x);
	int n = x->size;
	int iter = 0;
	double lambda;
	double fx;
	double fy;
	double uy;
    double fxs;
	double sy;
	double s_grad;
	const double dx = 1.0/524288;//sqrt(DBL_EPSILON);

	gsl_matrix* B = gsl_matrix_alloc(n,n);
	gsl_matrix* as = gsl_matrix_alloc(n,n);
	gsl_matrix* sa = gsl_matrix_alloc(n,n);
	gsl_matrix* dB = gsl_matrix_alloc(n,n);
	gsl_matrix_set_identity(B);
    gsl_vector* grad_F = gsl_vector_alloc(n);         
    gsl_vector* grad_Fs = gsl_vector_alloc(n);         
    gsl_vector* Dx = gsl_vector_alloc(n);    
    gsl_vector* y = gsl_vector_alloc(n);            
    gsl_vector* By = gsl_vector_alloc(n);    
    gsl_matrix* a = gsl_matrix_alloc(n,n);    
    gsl_matrix* s = gsl_matrix_alloc(n,n);    
    gsl_vector* u = gsl_vector_alloc(n);    
    gsl_matrix* R = gsl_matrix_alloc(n,n);
   
    while(1){
	//vector_print("x i while løkke qnewton",x);
	iter += 1;
	if (iter>2000) break;
        fx = f(x); 


//        for (int j =0; j<n; j++) {
//		gsl_vector_set(x,j,gsl_vector_get(x,j)+dx);
//		dfx = f(x);
//		df = dfx - fx;
//		gsl_vector_set(grad_F,j,df/dx);             
//		gsl_vector_set(x,j,gsl_vector_get(x,j)-dx);  
//        }

	numeric_gradient(f,x,grad_F);
	
	
	if(gsl_blas_dnrm2(grad_F)<eps){
		fprintf(stderr,"qnewton: |grad|<eps\n");
		printf("break1\n");
		break;
		}

	gsl_blas_dgemv(CblasNoTrans, -1, B, grad_F, 0, Dx);

    if(gsl_blas_dnrm2(Dx) < dx) {
		fprintf(stderr,"qnewton: |Dx|=%g < dx=%g\n",gsl_blas_dnrm2(Dx),dx);
		printf("break2\n");
        break;
        }

        lambda = 1;
	while(1){
		gsl_vector_memcpy(y,x);
		gsl_vector_add(y,Dx);
		fy = f(y);
		gsl_blas_ddot(Dx, grad_F, &s_grad);
		if(fy<fx+1e-4*s_grad){
			fprintf(stderr,"qnewton: fy=%g < fx=%g\n",fy,fx);
			break;
		}
		if (lambda < 1.0/32) {
			fprintf(stderr,"qnewton: lambda<dx !!!\n");
			gsl_matrix_set_identity(B);
			break;
		}
            lambda/=2;
            gsl_vector_scale(Dx,0.5);
        }

	// Opdatering af B til B + Î´B
//	gsl_vector_scale(Dx,lambda); // s-vektor bestemmes
	gsl_vector_add(Dx, x); 	// s <- s + x
	fxs = f(Dx); 			// funktionsvÃ¦rdien i (s+x)-vektoren


//	for (int j =0; j<n; j++) {
//		gsl_vector_set(Dx,j,gsl_vector_get(Dx,j)+dx);
//        dfx = f(Dx); 		// funktionsvÃ¦rdi i (s+x) + dx
//		df = dfx - fxs; 	// forskel i fkt.-vÃ¦rdi fÃ¸r og efter dx er lagt til
//        gsl_vector_set(grad_Fs,j,df/dx); //bestemmer gradienten af s+x            
//		gsl_vector_set(Dx,j,gsl_vector_get(Dx,j)-dx);   //gÃ¥r dx tilbage
//    	}
//
	numeric_gradient(f,Dx,grad_Fs);
	
	// u og y-vektor bestemmes
	gsl_vector_sub(Dx, x); // Vektoren Dx er lig s. (s <- s - x)	
	gsl_vector_memcpy(u,Dx); // u = s
	gsl_vector_sub(grad_Fs, grad_F); // y-vektor bestemmes y = grad_Fs
	gsl_blas_dgemv(CblasNoTrans, 1.0, B, grad_Fs, 0.0, By); // B*y -> By
	gsl_vector_sub(u, By); // u = s - By
	// prikprodukter uy og sy bestemmes
	gsl_blas_ddot(u, grad_Fs, &uy); // prikprodukt u^T*y = uy
	gsl_blas_ddot(Dx, grad_Fs, &sy); // prikprodukt s^T*y = sy
	
	if (fabs(sy)>eps){ // laver B om hvis sy er stÃ¸rre end eps
		gsl_blas_dger(1/sy,u,Dx,B);
// 
// gamma = uy/(2.0*sy); // gamma bestemmes 
// gsl_vector_scale(Dx, -gamma); // Dx <- -s*gamma
// gsl_vector_add(u, Dx);		  // u <- u + -s*gmma
// gsl_vector_scale(u, 1.0/sy);  // a bestemmes: u <- u/sy (a = u)
// 
// // laver matrixer til at bestemme Î´B
// for (int i = 0; i<u->size; i++){
// 	gsl_matrix_set(a,i,0,gsl_vector_get(u,i)); //a-matricen (n*1)
// 	gsl_matrix_set(s,i,0,gsl_vector_get(Dx,i)); //s-matricen (n*1)
// }
// gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, a, s, 0.0, as);
// gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, s, a, 0.0, sa);
// 
// gsl_matrix_memcpy(dB, as); 
// gsl_matrix_add(dB,sa); // as^T + sa^T = Î´B
// gsl_matrix_add(B,dB); // B + Î´B - B opdateres.
// 
	}
	
	// opdaterer til det nye step
	//vector_print("x",x);
	//vector_print("y",y);
	gsl_vector_memcpy(x,y); 
	//vector_print("x",x);
    fx = fy;
    }
    printf("iterations = %d\n", iter);
	
    gsl_matrix_free(B);
    gsl_matrix_free(as);
    gsl_matrix_free(sa);
    gsl_matrix_free(dB);
    gsl_vector_free(grad_F);
    gsl_vector_free(grad_Fs);
    gsl_vector_free(y);  
    gsl_vector_free(By); 
    gsl_vector_free(Dx);
    gsl_matrix_free(a);
    gsl_matrix_free(s);
    gsl_vector_free(u);
    gsl_matrix_free(R);
   
}

/*
#define TINY (1.0/524288)


void numeric_gradient
(double beta(gsl_vector*), gsl_vector*x, gsl_vector*grad){
double fx=beta(x);
for(int i=0;i<x->size;i++){
	double xi=gsl_vector_get(x,i);
	double dx=fabs(xi)*TINY;
	if(fabs(xi)<sqrt(TINY)) dx=TINY;
	gsl_vector_set(x,i,xi+dx);
	gsl_vector_set(grad,i,(beta(x)-fx)/dx);
	gsl_vector_set(x,i,xi);
	}
}

int qnewton(double beta(gsl_vector*), gsl_vector*x, double acc) {
int n=x->size,nsteps=0,nbad=0,ngood=0;
gsl_matrix* B=gsl_matrix_alloc(n,n);
gsl_vector* gx=gsl_vector_alloc(n);
gsl_vector* Dx=gsl_vector_alloc(n);
gsl_vector* z=gsl_vector_alloc(n);
gsl_vector* gz=gsl_vector_alloc(n);
gsl_vector* y=gsl_vector_alloc(n);
gsl_vector* u=gsl_vector_alloc(n);
gsl_matrix_set_identity(B);
numeric_gradient(beta,x,gx);
double fx=beta(x),fz;
while(nsteps<2000){
	nsteps++;
if(fx<acc){fprintf(stderr,"qnewton:converged: fx<acc=%g\n",acc); break;}
	gsl_blas_dgemv(CblasNoTrans,-1,B,gx,0,Dx);
//	if(gsl_blas_dnrm2(Dx)<TINY*gsl_blas_dnrm2(x))
//		{fprintf(stderr,"qnewton: |Dx|<TINY*|x|\n"); break;}
	if(gsl_blas_dnrm2(gx)<acc)
		{fprintf(stderr,"qnewton: |grad|<acc\n"); break;}
	double lambda=1;
	while(1){
		gsl_vector_memcpy(z,x);
		gsl_vector_add(z,Dx);
		fz=beta(z);
		double sTg; gsl_blas_ddot(Dx,gx,&sTg);
		if(fz<fx+0.01*sTg){
			ngood++;
			break;
			}
		if(lambda<1.0/32){
			nbad++;
			break;
			}
		lambda*=0.5;
		gsl_vector_scale(Dx,0.5);
		}
	numeric_gradient(beta,z,gz);
	gsl_vector_memcpy(y,gz);
	gsl_blas_daxpy(-1,gx,y); 
	gsl_vector_memcpy(u,Dx);
	gsl_blas_dgemv(CblasNoTrans,-1,B,y,1,u);
	double sTy,uTy;
	gsl_blas_ddot(Dx,y,&sTy);
	if(fabs(sTy)>1e-6){
		gsl_blas_ddot(u,y,&uTy);
		double gamma=uTy/2/sTy;
		gsl_blas_daxpy(-gamma,Dx,u);
		gsl_blas_dger(1.0/sTy,u,Dx,B);
		gsl_blas_dger(1.0/sTy,Dx,u,B);
		}
	gsl_vector_memcpy(x,z);
	gsl_vector_memcpy(gx,gz);
	fx=fz;
	}
gsl_matrix_free(B);
gsl_vector_free(gx);
gsl_vector_free(Dx);
gsl_vector_free(z);
gsl_vector_free(gz);
gsl_vector_free(y);
gsl_vector_free(u);
fprintf(stderr,"qnewton: nsteps=%i ngood=%i nbad=%i fx=%.1e\n"
		,nsteps,ngood,nbad,fx);
return nsteps;
}
*/
		