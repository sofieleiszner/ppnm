#include <gsl/gsl_integration.h>

double f (double t, void * params); 

double myerror(double z){
	gsl_function F;
	F.function=&f;
	F.params=(void*)&z;
	int limit=999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc (limit);
	double a=0, epsabs = 1e-6, epsrel = 1e-6, result, error;
    int  key = 6; 
    gsl_integration_qag(&F,a,z,epsabs,epsrel,limit,key,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}