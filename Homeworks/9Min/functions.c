#include <math.h>
#include <float.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>

static const double DELTA=sqrt(DBL_EPSILON);
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
	double f(gsl_vector* x), /* objective function */
	gsl_vector* x, /* on input: starting point, on exit: approximation to root */
	double eps /* accuracy goal, on exit |gradient| should be <eps */
	){
	int n = x->size;
	int iter = 0;
	double lambda;
	double fx;
	double fy;
	double uy;
	double sy;
	double s_grad;
	const double dx = sqrt(DBL_EPSILON);

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
	iter += 1;
	if (iter>1e6) break;
        fx = f(x); 

	numeric_gradient(f,x,grad_F);

	if(gsl_blas_dnrm2(grad_F)<eps){
		//fprintf(stderr,"qnewton: |grad|<eps\n");
		break;
		}

	gsl_blas_dgemv(CblasNoTrans, -1, B, grad_F, 0, Dx);

        if(gsl_blas_dnrm2(Dx) < dx) {
		//fprintf(stderr,"qnewton: |Dx|=%g < dx=%g\n",gsl_blas_dnrm2(Dx),dx);
                break;
                }

        lambda = 1;
	while(1){
		gsl_vector_memcpy(y,x);
		gsl_vector_add(y,Dx);
		fy = f(y);
		gsl_blas_ddot(Dx, grad_F, &s_grad);
		if(fy<fx+1e-4*s_grad){
			//fprintf(stderr,"qnewton: fy=%g < fx=%g\n",fy,fx);
			break;
		}
		if (lambda < dx) {
			//fprintf(stderr,"qnewton: lambda<dx !!!\n");
			gsl_matrix_set_identity(B);
			break;
		}
            lambda/=2;
            gsl_vector_scale(Dx,0.5);
        }

	// Opdatering af B til B + Î´B
	gsl_vector_add(Dx, x); 	// s <- s + x
	//fxs = f(Dx); 			// funktionsvÃ¦rdien i (s+x)-vektoren

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
	}
	
	// opdaterer til det nye step
	gsl_vector_memcpy(x,y); 
    fx = fy;
    }
    printf("Number of iterations = %d\n", iter);
	
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

void reflection(double* highest, double* centroid, int dim, double* reflected){
	for (int i=0; i<dim; i++) reflected[i]=2*centroid[i]-highest[i];
}

void expansion(double* highest, double* centroid, int dim, double* expanded){
	for (int i=0; i<dim; i++) expanded[i]=3*centroid[i]-2*highest[i];
}

void contraction(double* highest, double* centroid, int dim, double* contracted){
	for (int i=0; i<dim; i++) contracted[i]=0.5*centroid[i]-0.5*highest[i];
}

void reduction(double** simplex, int dim, int lo){
	for (int k=0; k<dim+1; k++) if (k!=lo) for (int i=0; i<dim; i++)
		simplex[k][i]=0.5*(simplex[k][i]+simplex[lo][i]);
}

double distance(double* a, double* b, int dim){
	double s=0; for (int i=0; i<dim; i++) s+=pow(b[i]-a[i], 2);
	return sqrt(s);
}

double size(double** simplex, int dim){
	double s=0; 
	for (int k=1; k<dim+1; k++){
		double dist = distance(simplex[0], simplex[k], dim);
		if (dist>s) s=dist; 
	}
	return s ;
}


void simplex_update(double** simplex , double* f_values, int d, int* hi, int* lo, double* centroid){
	* hi=0; * lo =0; 
	double highest=f_values[0], lowest=f_values[0];
	for (int k=1; k<d+1; k++){
		double next=f_values[k];
		if (next>highest){ 
			highest=next;  
			* hi=k; 
		}
		if (next<lowest){ 
			lowest=next ; 
			* lo=k; 
		} 
	}
	for (int i=0; i<d; i++){
		double s=0; 
		for (int k=0; k<d+1; k++) if (k!= *hi) s+=simplex[k][i];
		centroid[i]=s/d; 
	}
}

void simplex_initiate(
double fun(double *), double** simplex, double* f_values, int d,
int* hi , int* lo , double* centroid){
	for (int k=0; k<d+1; k++) f_values[k] = fun(simplex[k]);
	simplex_update(simplex, f_values, d, hi, lo, centroid);
}


//table 5
int downhill_simplex(double F(double *), double** simplex, int d, double simplex_size_goal){
int hi, lo, k=0; 
double centroid[d], F_value[d+1], p1[d], p2[d]; //opretter arrays - måske laves til vektorer?
simplex_initiate(F, simplex, F_value, d, &hi, &lo, centroid);
while (size(simplex, d) > simplex_size_goal){
	simplex_update(simplex, F_value, d, &hi, &lo, centroid);
	reflection(simplex[hi], centroid, d, p1); 
	double f_re=F(p1);
	if (f_re < F_value[lo]){ 	//reflection looks good : try expansion
		expansion(simplex[hi], centroid, d, p2); 
		double f_ex = F(p2);
		if (f_ex < f_re){ 		//accept expansion
			for (int i=0; i<d; i++){ simplex[hi][i]=p2[i]; F_value[hi] = f_ex;}
		}
		else { 					//reject expansion and accept reflection
			for (int i=0; i<d; i++){ simplex[hi][i]=p1[i]; F_value[hi]= f_re; }
		}
	}
	else { 						//reflection wasn't good
		if (f_re < F_value[hi]){//ok, accept reflection
			for (int i=0; i<d; i++){ simplex[hi][i]=p1[i]; F_value[hi]=f_re;}
		}
		else { 					//try contraction
			contraction(simplex[hi], centroid, d, p1); 
			double f_co = F(p1);
			if (f_co < F_value[hi]){ //accept contraction
				for (int i=0; i<d; i++){ simplex[hi][i]=p1[i]; F_value[hi]= f_co ;}
			}
			else { 				//do reduction
				reduction(simplex, d, lo);			
				simplex_initiate(F, simplex, F_value, d, &hi, &lo, centroid);
			}
		}
	}	
	k++;
}
printf("Minimum (x,y) = (%g, %g)\n", simplex[lo][0], simplex[lo][1]); 
return k;
}
			