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
    gsl_vector_set(dydx,0,gsl_vector_get(y,6));
    gsl_vector_set(dydx,1,gsl_vector_get(y,7));
    gsl_vector_set(dydx,2,gsl_vector_get(y,8));
    gsl_vector_set(dydx,3,gsl_vector_get(y,9));
    gsl_vector_set(dydx,4,gsl_vector_get(y,10));
    gsl_vector_set(dydx,5,gsl_vector_get(y,11));
    
    double y0 = gsl_vector_get(y,0);
    double y1 = gsl_vector_get(y,1);
    double y2 = gsl_vector_get(y,2);
    double y3 = gsl_vector_get(y,3);
    double y4 = gsl_vector_get(y,4);
    double y5 = gsl_vector_get(y,5);
    
    double R1 = pow(pow(y2-y0,2)+pow(y3-y1,2),0.5); //Længe^2 mellem r1 og r2
    double R2 = pow(pow(y4-y0,2)+pow(y5-y1,2),0.5); //Mellem r1 og r3
    double R3 = pow(pow(y4-y2,2)+pow(y5-y3,2),0.5); //r2 og r3
      
    
    gsl_vector_set(dydx,6,-(y0-y2)/pow(R1,3)-(y0-y4)/pow(R2,3)); 
    gsl_vector_set(dydx,7,-(y1-y3)/pow(R1,3)-(y1-y5)/pow(R2,3)); 
    gsl_vector_set(dydx,8, -(y2-y4)/pow(R3,3)-(y2-y0)/pow(R1,3)); 
    gsl_vector_set(dydx,9, -(y3-y5)/pow(R3,3)-(y3-y1)/pow(R1,3)); 
    gsl_vector_set(dydx,10,-(y4-y0)/pow(R2,3)-(y4-y2)/pow(R3,3)); 
    gsl_vector_set(dydx,11,-(y5-y1)/pow(R2,3)-(y5-y3)/pow(R3,3)); 

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
    n = 12; 
    gsl_vector* yaNewton = gsl_vector_alloc(n); 
    gsl_vector* ybNewton = gsl_vector_alloc(n); 
    gsl_vector_set(yaNewton,0, 0.97000436); 
    gsl_vector_set(yaNewton,1, -0.24308753); 
    gsl_vector_set(yaNewton,2, -0.97000436); 
    gsl_vector_set(yaNewton,3, 0.24308753); 
    gsl_vector_set(yaNewton,4, 0); 
    gsl_vector_set(yaNewton,5, 0); 
    gsl_vector_set(yaNewton,6, 0.93240737/2); 
    gsl_vector_set(yaNewton,7, 0.86473146/2); 
    gsl_vector_set(yaNewton,8, 0.93240737/2); 
    gsl_vector_set(yaNewton,9,0.86473146/2); 
    gsl_vector_set(yaNewton,10, -0.93240737); 
    gsl_vector_set(yaNewton,11, -0.86473146); 
    a = 0, b = 10, h = 0.01, acc = 1e-3, eps = 1e-3;
    driver(fNewton, a, yaNewton, b, ybNewton, h, acc, eps, "TheWayNewton.txt"); 
    
    
    
    return 0; 
} 