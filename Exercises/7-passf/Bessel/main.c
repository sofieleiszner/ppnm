#include<stdio.h>
#include<math.h>
#include <gsl/gsl_integration.h>

//Dette er måske ikke sådan 

double f (double t, void * params) {
  double n = 2;
  double x = *(double*)params;
  double f = 1/(M_PI)*cos(n*t-x*sin(t));
  return f;
}

double mybessel(double x){
	gsl_function F;
	F.function=&f;
	F.params=(void*)&x;
	int limit=999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc (limit);
	double a = 0, b = M_PI, epsabs = 1e-6, epsrel = 1e-6, result, error;
    int  key = 6; 
    gsl_integration_qag(&F,a,b,epsabs,epsrel,limit,key,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

int main(){
	for(double x=1.0;x<=32;x+=1.0/8)
		printf("%10g %10g\n",x,mybessel(x));
return 0;
}

