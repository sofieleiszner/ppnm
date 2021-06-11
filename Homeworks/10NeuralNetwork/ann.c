#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<math.h>
#include"ann.h"
double CP(gsl_vector * params);
double cost_function(gsl_vector* p);

void qnewton(
	double f(gsl_vector* x), /* objective function */
	gsl_vector* x, /* on input: starting point, on exit: approximation to root */
	double eps /* accuracy goal, on exit |gradient| should be <eps */
	);

ann*   ann_alloc(int n,double(*f)(double)){
	ann* network=malloc(sizeof(ann));
	gsl_vector * params=gsl_vector_alloc(n*3);
	network->n=n;
	network->f=f;
	network->params=params;
	return network;
}

void   ann_free(ann* network){
	gsl_vector_free(network->params);
	free(network);	
}

double ann_response(ann* network,double x){
	double sum=0;
	for(int i=0; i<network->n;i++){
		double ai=gsl_vector_get(network->params,i*3);
		double bi=gsl_vector_get(network->params,i*3+1);
		double wi=gsl_vector_get(network->params,i*3+2);
		sum+=network->f((x-ai)/bi)*wi;
		}
	return sum;
}



double ann_gradient(ann* network,double x){
	double sum=0;
	for(int i=0; i<network->n;i++){
		double ai=gsl_vector_get(network->params,i*3);
		double bi=gsl_vector_get(network->params,i*3+1);
		double wi=gsl_vector_get(network->params,i*3+2);
		sum+=wi*(exp(-(x-ai)*(x-ai)/(bi*bi))/bi-2.0*(x-ai)*(x-ai)*exp(-(x-ai)*(x-ai)/(bi*bi))/bi/bi/bi);
		}
	return sum;
}


double ann_integral(ann* network,double x){
	double sum=0;
	for(int i=0; i<network->n;i++){
		double ai=gsl_vector_get(network->params,i*3);
		double bi=gsl_vector_get(network->params,i*3+1);
		double wi=gsl_vector_get(network->params,i*3+2);
		//sum+=activationfunctionintegral((x-ai)/bi)*wi;
		sum+=-0.5*bi*exp(-(ai-x)*(ai-x)/(bi*bi))*wi;
		}
	return sum;
}


void ann_train(ann* network,gsl_vector* xfunc,gsl_vector* yexact){
	double eps=0.01;
	gsl_vector *p=gsl_vector_alloc(network->params->size);
	gsl_vector_memcpy(p,network->params);
	qnewton(cost_function,p,eps);
	gsl_vector_memcpy(network->params,p);
	gsl_vector_free(p);

}

/*
void ann_trainC(ann* network,gsl_vector* xfunc,gsl_vector* yexact){
	double eps=0.01;
	gsl_vector *p=gsl_vector_alloc(network->params->size);
	gsl_vector_memcpy(p,network->params);
	qnewton(cost_function,p,eps);
	gsl_vector_memcpy(network->params,p);
	gsl_vector_free(p);
}
*/

