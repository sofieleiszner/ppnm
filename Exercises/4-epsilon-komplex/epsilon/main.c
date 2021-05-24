#include <stdio.h>
#include <limits.h>
#include <float.h>
 


int equall(double a, double b, double tau, double epsilon);


int main(){
    
    printf("Task 1\n");
    
    printf("Task 1i\n");
	
    int i=1; 
    while(i+1>i) {i++;}
	printf("My max int (while loop) = %i \n", i);
    
    i=1;
	do{
	i++;}
	while(i+1>i);
	printf("My max int (do while loop) = %i\n", i);
    
    for(i=1; i<i+1; i++){}
	printf("My max int (for loop) = %i\n", i);
    printf("Comparisson with INT_MAX from limits.h, INT_MAX = %i\n", INT_MAX); 
  
    
    printf("\nTask 1ii\n");
	i=1; 
    while(i-1<i) {i--;}
	printf("My minimum int (while loop) = %i \n", i);
    
    i=1;
	do{
	i--;}
	while(i-1<i);
    printf("My minimum int (do while loop) = %i\n", i);
    
    for(i=1; i-1<i; i--){}
    printf("My minimum int (for loop) = %i\n", i);

    printf("Comparisson with INT_MIN from limits.h, INT_MIN = %i\n", INT_MIN); 
    
    
    printf("\nTask 1iii\n"); 
    printf("With while loops:\n"); 
    float y  = 1; 
    while(1+y != 1) {y/=2;} y*=2; 
	printf("With float, machine eps = %.10f\n", y); 
	double x = 1; 
    while(1+x != 1) {x/=2;} x*=2;
	printf("With double, machine eps = %.10g\n", x);
	long double z = 1; 
    while(1+z != 1) {z/=2;} 
    z*=2;
	printf("With long double, machine eps = %.10Lg\n", z);  
    
    printf("With do while loops:\n"); 
    y  = 1; 
    do{y/=2;}
    while (1+y != 1);
    y*=2; 
	printf("With float, machine eps = %.10f\n", y); 
	x = 1; 
        do{x/=2;}
    while (1+x != 1);
    x*=2; 
	printf("With double, machine eps = %.10g\n", x);
	z = 1; 
    do{z/=2;}
    while (1+z != 1);
    z*=2; 
    printf("With long double, machine eps = %.10Lg\n", z);  
    
    printf("With for loops:\n"); 
    for(y=1; 1+y != 1; y/=2){} 
    y*=2; 
	printf("With float, machine eps = %.10f\n", y); 
    for(x=1; 1+x != 1; x/=2){}
    x*=2;
	printf("With double, machine eps = %.10g\n", x);
    for(z=1; 1+z != 1; z/=2){}
    z*=2;
	printf("With long double, machine eps = %.10Lg\n", z);  

    printf("Comparisson with FLT_EPSILON from limits.h,  FLT_EPSILON = %.10f\n", FLT_EPSILON); 
    printf("Comparisson with DBL_EPSILON from limits.h,  DBL_EPSILON = %.10g\n", DBL_EPSILON); 
    printf("Comparisson with LDBL_EPSILON from limits.h, LDBL_EPSILON = %.10Lg\n", LDBL_EPSILON); 
    
    
    
    printf("\nTask 2\n");
    printf("\nTask 2i\n");

    int max = 2147483647/3; 
	float sum_up_float = 0; 
	i = 1;
	for (i = 1; i<max; i++) {
		sum_up_float += 1.0f/i;
	}
	printf("sum_up_float = %g \n", sum_up_float);
	float sum_down_float = 1.0f/max;
	for (i = 1; i<max; i++) {
	sum_down_float +=1.0f/(max-i);
	}
	printf("sum_down_float = %g \n", sum_down_float); 
    
    printf("It gives greater precision to add the small numbers first.  \n\
If you first add large numbers and then add the small numbers, the small \n\
numbers will just be some insignificant small digits that cannot be saved due \n\
to limited number of bits (in computer memory). \n");

    printf("\nTask 2iii\n");
    printf("No, the harmonic series is divergent\n"); 

    printf("\nTask 2iv\n");
    double sum_up_double = 0;
    for (i = 1; i<max; i++) {
		sum_up_double += 1.0f/i;
	}
	printf("sum_up_double = %.10g \n", sum_up_double);
	double sum_down_double = 1.0f/max;
	for (i = 1; i<max; i++) {
	sum_down_double +=1.0f/(max-i);
	}
	printf("sum_down_double = %.10g \n", sum_down_double); 
    printf("The difference between the two sums is now very small, because of the increased precision of doubles compared to floats\n"); 
    
    printf("\nTask 3\n");

    printf("Checks if a and b are equal with regards to the relative precision (eps) and \
absolute precision (tau) given to the function \n"); 
    printf("The function returns 1 if the numbers 'a' and 'b' are equal with absolute precision 'tau',\
or are equal with relative precision 'epsilon', and returns 0 otherwise.\n");
    
    
    int equal_result = equall(8.3, 8.2, 0.5, 0.1);
	printf("equal(a = 8.3, b = 8.2, tau = 0.5, eps = 0.1) = %d \n", equal_result);
    equal_result = equall(8.3, 9.0, 0.5, 0.01);
    printf("equal(a = 8.3, b = 9.0, tau = 0.5, eps = 0.01) = %d \n", equal_result);
return 0;
}
