#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>




void rkstep12(
	void f(double x,gsl_vector* y,gsl_vector* dydx), /* the f from dy/dt=f(t,y) */
	double x,              /* the current value of the variable */
	gsl_vector* yx,            /* the current value y(t) of the sought function */
	double h,              /* the step to be taken */
	gsl_vector* yh,             /* output: y(t+h) */
	gsl_vector* dy            /* output: error estimate */
    ); 
    
    
void driver(
	void f(double x,gsl_vector* y ,gsl_vector* dydx), /* right-hand-side of dy/dt=f(t,y) */
	double a,                     /* the start-point a */
	gsl_vector* ya,                     /* y(a) */
	double b,                     /* the end-point of the integration */
	gsl_vector* yb,                     /* y(b) to be calculated */
	double h,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps,                    /* relative accuracy goal */
    char filename[]
);


void f1(double x, gsl_vector* y, gsl_vector* dydx){
    gsl_vector_set(dydx,0,gsl_vector_get(y,1)); 
    gsl_vector_set(dydx,1,-gsl_vector_get(y,0));
    
}

void fSIR(double x, gsl_vector* y, gsl_vector* dydx){
    int N = 1000; 
    int Tc = 4; 
    int Tr = 10; 
    
    gsl_vector_set(dydx,0,-gsl_vector_get(y,1)*gsl_vector_get(y,0)/(N*Tc));
    gsl_vector_set(dydx,1,gsl_vector_get(y,1)*gsl_vector_get(y,0)/(N*Tc)-gsl_vector_get(y,1)/Tr);
    gsl_vector_set(dydx,2,gsl_vector_get(y,1)/Tr);
    
}


void fNewton(double x, gsl_vector* y, gsl_vector* dydx){
    double G = 6.67*pow(10,-11); //N*m^2/kg^2 
    double m1 = 1e12; //kg
    double m2 = 2e12; //kg
    double m3 = 3e12; //kg
    double r1 = gsl_vector_get(y,0);
    double r2 = gsl_vector_get(y,1);
    double r3 = gsl_vector_get(y,2);
    gsl_vector_set(dydx,0,-G*m2*(r1-r2)/pow((r1-r2),3)-G*m3*(r1-r3)/pow((r1-r3),3)); 
    gsl_vector_set(dydx,1,-G*m3*(r2-r3)/pow((r2-r3),3)-G*m1*(r2-r1)/pow((r2-r1),3)); 
    gsl_vector_set(dydx,2,-G*m1*(r3-r1)/pow((r3-r1),3)-G*m2*(r3-r2)/pow((r3-r2),3)); 
}


int main() {
    //Checks sinusfunction
    int n = 2; 
    gsl_vector* ya = gsl_vector_alloc(n); 
    gsl_vector* yb = gsl_vector_alloc(n); 
    gsl_vector_set(ya,0,0); 
    gsl_vector_set(ya,1,1);   
    double a = 0, b = 7, h = 0.1, acc = 1e-2, eps = 1e-2;
    driver(f1, a, ya, b, yb, h, acc, eps, "TheWay.txt"); 
    
    //SIR
    n = 3; 
    gsl_vector* yaSIR = gsl_vector_alloc(n); 
    gsl_vector* ybSIR = gsl_vector_alloc(n); 
    gsl_vector_set(yaSIR,0,900); 
    gsl_vector_set(yaSIR,1,10); 
    gsl_vector_set(yaSIR,2,0); 
    a = 0, b = 150, h = 0.1, acc = 1e-3, eps = 1e-3;
    driver(fSIR, a, yaSIR, b, ybSIR, h, acc, eps, "TheWaySIR.txt"); 
    
    //Newton
    n = 3; 
    gsl_vector* yaNewton = gsl_vector_alloc(n); 
    gsl_vector* ybNewton = gsl_vector_alloc(n); 
    gsl_vector_set(yaNewton,0,1); 
    gsl_vector_set(yaNewton,1,2); 
    gsl_vector_set(yaNewton,2,4); 
    a = 0, b = 100, h = 0.01, acc = 1e-3, eps = 1e-3;
    driver(fNewton, a, yaNewton, b, ybNewton, h, acc, eps, "TheWayNewton.txt"); 
    
    
    
    return 0; 
} 