#include<gsl/gsl_vector.h>

int binsearch(int n,gsl_vector *x,double Z);
typedef struct {gsl_vector *x, *y, *b, *c;} qspline;

qspline* qspline_alloc(gsl_vector *x, gsl_vector* y){
	qspline *s=(qspline*)malloc(sizeof(qspline));
	s->b=gsl_vector_alloc((x->size)-1);
	s->c=gsl_vector_alloc((x->size)-1);
	s->x=gsl_vector_alloc(x->size);
	s->y=gsl_vector_alloc(x->size);
	s->x=x; s->y=y;
	int i; gsl_vector * h=gsl_vector_alloc((x->size)-1); gsl_vector * p=gsl_vector_alloc((x->size)-1);

	for(i=0;i<(x->size)-1;i++){
	gsl_vector_set(h,i,gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
	gsl_vector_set(p,i,(gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/gsl_vector_get(h,i));
	}
	gsl_vector_set(s->c,0,0);

	for(i=0;i<(x->size)-2;i++){
	gsl_vector_set(s->c,i+1,(gsl_vector_get(p,i+1)-gsl_vector_get(p,i)-gsl_vector_get(s->c,i)*gsl_vector_get(h,i))/gsl_vector_get(h,i+1));
	}
	gsl_vector_set(s->c,x->size-2,gsl_vector_get(s->c,i-2)/2);

	for(i=x->size-3;i>=0;i--){gsl_vector_set(s->c,i,(gsl_vector_get(p,i+1)-gsl_vector_get(p,i)-gsl_vector_get(s->c,i+1)*gsl_vector_get(h,i+1))/gsl_vector_get(h,i));
	}
	for(i=0;i<x->size-1;i++){gsl_vector_set(s->b,i,gsl_vector_get(p,i)-gsl_vector_get(s->c,i)*gsl_vector_get(h,i));
	}
return s;
}

double qspline_eval(qspline *s, double z){
	int index=binsearch(s->x->size,s->x,z);
	double h=z-gsl_vector_get(s->x,index);
	return gsl_vector_get(s->y,index)+h*(gsl_vector_get(s->b,index)+h*gsl_vector_get(s->c,index));
}
void quad_spline_free(qspline *s){
	gsl_vector_free(s->x); gsl_vector_free(s->y); gsl_vector_free(s->b); gsl_vector_free(s->x);
}
