#include <stdio.h> 
#include <math.h> 
#include <complex.h>


int main(){
	
    printf("Subtask 1:\n");
    int n = 5;
	double x = 0.5; 
	complex z = csqrt(-2); 
	complex im = csqrt(-1);
	printf("Gamma function, Γ(5) =%g\n", gamma(n));
	printf("Bessel function, J_1(0.5) = %g \n", j1(x)); 
	printf("sqrt(2) = %g I%g \n", creal(z), cimag(z)); 
	printf("exp(i*pi) = %g I%g \n", creal(cexp(im*M_PI)), cimag(cexp(im*M_PI))); 
	printf("exp(i) = %g I%g \n", creal(cexp(im)), cimag(cexp(im))); 
	printf("i^e = %g I%g \n", creal(cpow(im,M_E)), cimag(cpow(im,M_E))); 
	printf("i^i = %g I%g \n", creal(cpow(im,im)), cimag(cpow(im,im))); 
    printf("\n");
    
    printf("Subtask 2:\n");
	float x_float = 1.f/9;
	double x_double = 1./9;
	long double x_long_double = 1.L/9;
	printf("Float = %.25g \n", x_float);
	printf("Double = %.25lg \n", x_double);
	printf("Long double = %.25Lg \n", x_long_double);
    printf("Float holds 9 significant digits\n");
    printf("Double holds 16 significant digits\n");
    printf("Long double holds 19 significant digits\n");


return 0;
} 

