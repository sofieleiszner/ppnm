#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_integration.h>

double infintegrate(double f(double x), double a, double b, double acc, double eps);

int calls=0;

double f1(double x){
			calls++;
			return sqrt(x);
		}
        
double f2(double x){
			calls++;
			return 4*sqrt(1-x*x);
		}        

double f3(double x){
			calls++;
			return 1/sqrt(x);
		}

double g(double x){return f3(cos(x))*sin(x);}; 	


double g_gsl(double x, void * params){
	return f3(cos(x))*sin(x);
}; 	

double f4(double x){calls++; return 2.0/(1.0 + x*x);};

double gsl_f4(double x, void * params){calls++; return 2.0/(1.0 + x*x);};
//double f2(double x, void * params){
//	double f2 = 4*sqrt(1-pow(x,2));
//	return f2;
//};





int main(){
	printf("Task A\n"); 
    
	double a=0.0, b=1.0; 
	double acc=0.001, eps=0.001;
    double exact = 2.0/3; 
	double Q1 = infintegrate(f1,a,b,acc,eps); 
    printf("Integration function infitegrate is implemented in functions.c\n\n"); 
	printf("Integration of sqrt(x) from %g to %g with eps = %g and accuracy = %g\n", a, b, eps, acc);
    printf("Calculated, Q = %.5g \n", Q1);
    printf("Exact = %g\n", exact); 
    printf("Number of call = %d\n", calls); 
    printf("Estimated error = %g\n", acc+fabs(Q1)*eps); 
    printf("Actual error = %g\n", fabs(Q1-exact)); 
    printf("\n"); 
    
    calls = 0; 
    exact = M_PI; 
	double Q1b = infintegrate(f2,a,b,acc,eps); 
	printf("Integration of 4*sqrt(1-x*x) from %g to %g with eps = %g and accuracy = %g\n", a, b, eps, acc);
    printf("Calculated, Q = %.5g \n", Q1b);
    printf("Exact = %g\n", exact); 
    printf("Number of call = %d\n", calls); 
    printf("Estimated error = %g\n", acc+fabs(Q1b)*eps); 
    printf("Actual error = %g\n", fabs(Q1b-exact)); 
    printf("\n"); 
    
    calls = 0; 
    exact = 2; 
	double Q1c = infintegrate(f3,a,b,acc,eps); 
	printf("Integration of 1/sqrt(x) from %g to %g with eps = %g and accuracy = %g\n", a, b, eps, acc);
    printf("Calculated, Q = %.5g \n", Q1c);
    printf("Exact = %g\n", exact); 
    printf("Number of call = %d\n", calls); 
    printf("Estimated error = %g\n", acc+fabs(Q1c)*eps); 
    printf("Actual error = %g\n", fabs(Q1c-exact)); 
    
	
	// Opgave B
    printf("\n---------------------------------------------\n");
	printf("Task B\n"); 
    printf("Implements Clenshaw Curtis variable transformation\n\n"); 
	calls = 0;
	a = 0; b = 1;
	printf("Integration of 1/sqrt(x) from %g to %g with eps = %g and accuracy = %g using CC variable transformation\n", a, b, eps, acc);
	a = acos(b); // Laver variable om
	b = acos(a); // Bytter om på grænserne fordi dx = -sinθdθ (se billede + afl-beskrivelse) - gælder uanset hvad a og b er (ikke inf)   
    printf("Calculates new limits: a = %g, b = %g\n", a, b);
	double Q2 = infintegrate(g,a,b,acc,eps);
    exact = 2;
    printf("Calculated, Q = %.5g \n", Q2);
    printf("Exact = %g\n", exact); 
    printf("Number of call = %d\n", calls); 
    printf("Estimated error = %g\n", acc+fabs(Q2)*eps); 
    printf("Actual error = %g\n", fabs(Q2-exact)); 
    printf("Clear improvement in number of calls\n\n"); 
    
    printf("Comparisson with GSL:\n"); 
	size_t size = 10000;
	double result, abserr;
	calls = 0;
	gsl_function F;
    F.function=&g_gsl;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(size);
	gsl_integration_qags(&F, a, b, acc, eps, size, w, &result, &abserr);
	printf("Calculated, Q = %.5g \n", result);
    printf("Exact = %g\n", exact); 
    printf("Number of call = %d\n", calls); 
    printf("Estimated error = %g\n", acc+fabs(result)*eps); 
    printf("Actual error = %g\n", fabs(result-exact)); 
    	
	//Opgave C
    printf("---------------------------------------------\n");
	printf("\nTask C: Infinite limits\n\n"); 
    calls = 0;
    a = 0; 
	b = INFINITY;
    exact = M_PI; 
    printf("Integration of 2.0/(1.0 + x*x) from %g to %g with eps = %g and accuracy = %g\n", a, b, eps, acc);
	double Q_inf = infintegrate(f4,a,b,acc,eps);
    printf("Calculated, Q = %.5g \n", Q_inf);
    printf("Exact = %g\n", exact); 
    printf("Number of call = %d\n", calls); 
    printf("Estimated error = %g\n", acc+fabs(Q_inf)*eps); 
    printf("Actual error = %g\n", fabs(Q_inf-exact)); 
    printf("\n");
    
    printf("Comparisson with GSL routine\n"); 
	calls = 0;
    F.function=&gsl_f4;
	gsl_integration_qagiu(&F, a, acc, eps, size, w, &result, &abserr);
	printf("Calculated, Q = %.5g \n", result);
    printf("Exact = %g\n", exact); 
    printf("Number of call = %d\n", calls); 
    printf("Estimated error = %g\n", acc+fabs(result)*eps); 
    printf("Actual error = %g\n", fabs(result-exact)); 	
return 0;
}



