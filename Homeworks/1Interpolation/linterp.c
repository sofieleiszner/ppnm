#include<math.h>
#include<gsl/gsl_vector.h>

int binsearch(int n, gsl_vector *x, double z);

double linterp(gsl_vector *x,gsl_vector *y, double Z){

	int interval=binsearch(x->size,x,Z);
	double xi=gsl_vector_get(x,interval);
	double xip1=gsl_vector_get(x,interval+1);
	double yi=gsl_vector_get(y,interval);
	double yip1=gsl_vector_get(y,interval+1);
	
	double yiz=yi+(yip1-yi)/(xip1-xi)*(Z-xi);
return yiz;
}
