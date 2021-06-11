#include<gsl/gsl_vector.h>
#include<math.h>

int binsearch(int n, gsl_vector *x, double z);
typedef struct {gsl_vector *x, *y, *b, *c;} qspline;

double quad_diff(qspline *s, double Z){
	int maxinterval=binsearch(s->x->size,s->x,Z);

	double xi=gsl_vector_get(s->x,maxinterval);
	double ci=gsl_vector_get(s->c,maxinterval);
	double bi=gsl_vector_get(s->b,maxinterval);
	double diff=bi+2*ci*(Z-xi);
return diff;
}
