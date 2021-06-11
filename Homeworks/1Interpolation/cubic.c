#include<assert.h>
#include<gsl/gsl_vector.h>
#include<math.h>

int binsearch(int n,gsl_vector *x,double Z);

typedef struct {gsl_vector *x, *y, *b, *c, *d;} cubic_spline;

cubic_spline* cubic_spline_alloc(gsl_vector *x, gsl_vector* y){
	cubic_spline* s=(cubic_spline*)malloc(sizeof(cubic_spline));
	int n = x->size;

	s->b=gsl_vector_alloc(n);
	s->c=gsl_vector_alloc(n-1);
	s->d=gsl_vector_alloc(n-1);
	s->x=gsl_vector_alloc(n);
	s->y=gsl_vector_alloc(n);
	s->x=x; s->y=y;

	gsl_vector * h=gsl_vector_alloc(n-1);
	gsl_vector * p=gsl_vector_alloc(n-1);

	for(int i =0; i < n-1; i++)
		{
		gsl_vector_set(h,i,gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
		assert(gsl_vector_get(h,i)>0);
		}

	for(int i =0; i < n-1; i++)
		{
		gsl_vector_set(p,i,(gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/gsl_vector_get(h,i));
		}

	gsl_vector * D=gsl_vector_alloc(n);
	gsl_vector * Q=gsl_vector_alloc(n-1);
	gsl_vector * B=gsl_vector_alloc(n);

	gsl_vector_set(D,0,2);
	for(int i=0; i < n-2; i++)
		{
		gsl_vector_set(D,i+1,2*gsl_vector_get(h,i)/gsl_vector_get(h,i+1)+2);
		}

	gsl_vector_set(D,n-1,2);
	gsl_vector_set(Q,0,1);

	for(int i=0; i < n-2; i++)
		{
		gsl_vector_set(Q,i+1,gsl_vector_get(h,i)/gsl_vector_get(h,i+1));
		}

	for(int i = 0; i < n - 2; i++)
		{
		gsl_vector_set(B, i+1, 3*(gsl_vector_get(p,i)+gsl_vector_get(p,i+1)*gsl_vector_get(h,i)/gsl_vector_get(h,i+1)));
		}

	gsl_vector_set(B,0,3*gsl_vector_get(p,0));
	gsl_vector_set(B, n-1,3*gsl_vector_get(p,n-2));

	for (int i =1; i<n; i++)
		{
		gsl_vector_set(D,i,gsl_vector_get(D,i)-gsl_vector_get(Q,i-1)/gsl_vector_get(D,i-1));
		gsl_vector_set(B,i,gsl_vector_get(B,i)-gsl_vector_get(B,i-1)/gsl_vector_get(D,i-1));
		}

	gsl_vector_set(s->b, n-1, gsl_vector_get(B, n-1)/gsl_vector_get(D,n-1));

	for(int i=n-2; i>=0; i--)
		{
		gsl_vector_set(s->b, i, (gsl_vector_get(B,i)-gsl_vector_get(Q,i)*gsl_vector_get(s->b,i+1))/gsl_vector_get(D,i));
		}
	for(int i =0; i<n-1; i++)
		{
		gsl_vector_set(s->c, i, (-2*gsl_vector_get(s->b,i)-gsl_vector_get(s->b,i+1)+3*gsl_vector_get(p,i))/gsl_vector_get(h,i));
		gsl_vector_set(s->d, i, (gsl_vector_get(s->b,i)+gsl_vector_get(s->b,i+1)-2*gsl_vector_get(p,i))/gsl_vector_get(h,i)/gsl_vector_get(h,i));
		}

return s;
}

double cubic_spline_eval(cubic_spline*s, double z) {
	int index = binsearch(s->x->size, s->x, z);
	double h = z - gsl_vector_get(s->x,index);
	return gsl_vector_get(s->y,index)+h*(gsl_vector_get(s->b,index)+h*(gsl_vector_get(s->c,index)+h*(gsl_vector_get(s->d,index))));
}

void cubic_spline_free(cubic_spline *s){
	gsl_vector_free(s->x); gsl_vector_free(s->y); gsl_vector_free(s->b); gsl_vector_free(s->c); gsl_vector_free(s->d); free(s);
}



double cubic_diff(cubic_spline*s, double Z){
	int maxinterval=binsearch(s->x->size,s->x,Z);
	double xi=gsl_vector_get(s->x,maxinterval);
	double ci=gsl_vector_get(s->c,maxinterval);
	double bi=gsl_vector_get(s->b,maxinterval);
	double di=gsl_vector_get(s->d,maxinterval);
	double diff=bi+2*ci*(Z-xi)+3*di*(Z-xi)*(Z-xi);
return diff;
}





double cubic_integ(cubic_spline *s, double Z){
	int maxinterval=binsearch(s->x->size,s->x,Z);
	double integral=0;

	for(int i=0;(i+1)<=maxinterval;i++)
		{
		double xi=gsl_vector_get(s->x,i);
		double xip1=gsl_vector_get(s->x,i+1);
		double yi=gsl_vector_get(s->y,i);
		double ai=yi;
		double bi=gsl_vector_get(s->b,i);
		double ci=gsl_vector_get(s->c,i);
		double di=gsl_vector_get(s->d,i);
		integral+=(double)ai*(xip1-xi)+(double)bi/2*pow((xip1-xi),2)+(double)ci/3*pow((xip1-xi),3)+(double)di/4*pow(xip1-xi,4);
		}

	double xi=gsl_vector_get(s->x,maxinterval);
	double yi=gsl_vector_get(s->y,maxinterval);
	double ai=yi;
	double bi=gsl_vector_get(s->b,maxinterval);
	double ci=gsl_vector_get(s->c,maxinterval);
	double di=gsl_vector_get(s->d,maxinterval);
	integral+=(double)ai*(Z-xi)+(double)bi/2*pow((Z-xi),2)+(double)ci/3*pow((Z-xi),3)+(double)di/4*pow(Z-xi,4);

return integral;
}