#include<stdio.h>
#include<math.h>


int equall(double a, double b, double tau, double epsilon) {
	if (fabs(a-b) < tau) {
        return 1;}
        
	if (fabs(a-b)/(fabs(a)+fabs(b)) < epsilon/2 ) {
        return 1;
        }
	return 0;}
