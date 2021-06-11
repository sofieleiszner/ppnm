#include<gsl/gsl_vector.h>
#include<math.h>

int binsearch(int n, gsl_vector *x, double z);

double linterp_integ(gsl_vector *x, gsl_vector *y, double Z){
	int maxinterval=binsearch(x->size,x,Z);
	double integral=0;
	for(int i=0;(i+1)<=maxinterval;i++)
		{
		double xi=gsl_vector_get(x,i);
		double xip1=gsl_vector_get(x,i+1);
		double yi=gsl_vector_get(y,i);
		double yip1=gsl_vector_get(y,i+1);
		double a=(yip1-yi)/(xip1-xi);
		double b=yi;
		integral+=(double)a/2*pow((xip1-xi),2)+b*(xip1-xi);
		}
	double xi=gsl_vector_get(x,maxinterval);
	double xip1=gsl_vector_get(x,maxinterval+1);
	double yi=gsl_vector_get(y,maxinterval);
	double yip1=gsl_vector_get(y,maxinterval+1);
	double a=(yip1-yi)/(xip1-xi);
	double b=yi;
	integral+=(double)a/2*pow((Z-xi),2)+b*(Z-xi);
return integral;
}
