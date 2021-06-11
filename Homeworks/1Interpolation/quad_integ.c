#include<gsl/gsl_vector.h>
#include<math.h>

typedef struct {gsl_vector *x, *y, *b, *c;} qspline;

int binsearch(int n, gsl_vector *x, double z);

double quad_integ(qspline *s, double Z){
	int maxinterval=binsearch(s->x->size,s->x,Z);
	double integral=0;

	for(int i=0;(i+1)<=maxinterval;i++)
		{
		double xi=gsl_vector_get(s->x,i);
		double xip1=gsl_vector_get(s->x,i+1);
		double yi=gsl_vector_get(s->y,i);
//		double yip1=gsl_vector_get(y,i+1);
		double ai=yi;
		double bi=gsl_vector_get(s->b,i);
		double ci=gsl_vector_get(s->c,i);
		integral+=(double)ai*(xip1-xi)+(double)bi/2*pow((xip1-xi),2)+(double)ci/3*pow((xip1-xi),3);
		}

	double xi=gsl_vector_get(s->x,maxinterval);
	double yi=gsl_vector_get(s->y,maxinterval);
	double ai=yi;
	double bi=gsl_vector_get(s->b,maxinterval);
	double ci=gsl_vector_get(s->c,maxinterval);
	integral+=(double)ai*(Z-xi)+(double)bi/2*pow((Z-xi),2)+(double)ci/3*pow((Z-xi),3);

return integral;
}

