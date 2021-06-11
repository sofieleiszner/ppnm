#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>

void plainmcA(double f(gsl_vector* x),gsl_vector* a,gsl_vector* b,int N,gsl_vector* results);
void plainmcB(double f(gsl_vector* x),gsl_vector* a,gsl_vector* b,int N,gsl_vector* results);

double f1(gsl_vector * x){
	return sqrt(gsl_vector_get(x,0));
	}

double f2(gsl_vector *x){
	if (1-cos(gsl_vector_get(x,0))*cos(gsl_vector_get(x,1))*cos(gsl_vector_get(x,2))==0) return 10;
	return 1/pow(M_PI,3)*pow(1-cos(gsl_vector_get(x,0))*cos(gsl_vector_get(x,1))*cos(gsl_vector_get(x,2)),-1);
	}

int main(){
	//A
    printf("Task A and B\n\n"); 
    printf("Plain MC and quasi-random MC integrators are implemented in functions.c\n");
    printf("Calculates integral of sqrt(x) from 0 to 1\n");
	int dim1=1;
	gsl_vector* a1=gsl_vector_alloc(dim1);
	gsl_vector* b1=gsl_vector_alloc(dim1);
	gsl_vector_set(a1,0,0);
	gsl_vector_set(b1,0,1);
	int N1=100000;
	gsl_vector* results1=gsl_vector_alloc(2);
    plainmcA(f1,a1,b1,N1,results1);
    printf("Integral of sqrt(x) from 0 to 1 with plain MC        = %.10g, error = %.10g\n",gsl_vector_get(results1,0),gsl_vector_get(results1,1));
    plainmcB(f1,a1,b1,N1,results1);
    printf("Integral of sqrt(x) from 0 to 1 with quasi-random MC = %.10g, error = %.10g\n",gsl_vector_get(results1,0),gsl_vector_get(results1,1));
    printf("Exact result: 2/3=0.66667 \n"); 
	FILE * file1=fopen("A.txt","w");
	FILE * file2=fopen("B.txt","w");
	for(int N=1000; N<300000; N+=5000){
		plainmcA(f1,a1,b1,N,results1);
		fprintf(file1,"%d %.25g %.25g\n",N,gsl_vector_get(results1,0),gsl_vector_get(results1,1));
		plainmcB(f1,a1,b1,N,results1);
		fprintf(file2,"%d %.25g %.25g\n",N,gsl_vector_get(results1,0),gsl_vector_get(results1,1));
		}
	fclose(file1);
	fclose(file2);
	
    printf("\nCalculates the long integral in the task description of task A.\n");
	int dim2=3;
	gsl_vector*a2=gsl_vector_alloc(dim2);
	gsl_vector*b2=gsl_vector_alloc(dim2);
	gsl_vector_set(a2,0,0);
	gsl_vector_set(b2,0,M_PI);
	gsl_vector_set(a2,1,0);
	gsl_vector_set(b2,1,M_PI);
	gsl_vector_set(a2,2,0);
	gsl_vector_set(b2,2,M_PI);
	int N2=100000;
	gsl_vector* results2=gsl_vector_alloc(2);
	plainmcA(f2,a2,b2,N2,results2);
	printf("Integral of long integral in task A with plain MC        = %.10g, error = %.10g\n",gsl_vector_get(results2,0),gsl_vector_get(results2,1));
	plainmcB(f2,a2,b2,N2,results2);
    printf("Integral of long integral in task A with quasi-random MC = %.10g, error = %.10g\n",gsl_vector_get(results2,0),gsl_vector_get(results2,1));
    printf("Exact result: 1.3932039296856768591842462603255\n"); 
	FILE * file3=fopen("AA.txt","w");
	FILE * file4=fopen("BB.txt","w");
	for(int N=1000; N<300000; N+=2500){
		plainmcA(f2,a2,b2,N,results2);
		fprintf(file3,"%d %.25g %.25g\n",N,gsl_vector_get(results2,0),gsl_vector_get(results2,1));
		
		plainmcB(f2,a2,b2,N,results2);
		fprintf(file4,"%d %.25g %.25g\n",N,gsl_vector_get(results2,0),gsl_vector_get(results2,1));
		}
	fclose(file3);
	fclose(file4);
	printf("\nFor comparison of errors, see the two plots: Error_Integral1.png and Error_Integral2.png\n");
    printf("The errors fall of much faster for the quasi-random MC than plain MC\n"); 
    
	//B
	
	//int dim3=3;
	//gsl_vector*a3=gsl_vector_alloc(dim3);
	//gsl_vector*b3=gsl_vector_alloc(dim3);
	//gsl_vector_set(a3,0,0);
	//gsl_vector_set(b3,0,M_PI);
	//gsl_vector_set(a3,1,0);
	//gsl_vector_set(b3,1,M_PI);
	//gsl_vector_set(a3,2,0);
	//gsl_vector_set(b3,2,M_PI);
	//int N3=1000;
	//gsl_vector* results3=gsl_vector_alloc(2);
	//plainmcB(f2,a3,b3,N3,results3);
	//printf("B integral=%.25g, error=%.25g\n",gsl_vector_get(results3,0),gsl_vector_get(results3,1));
	
	//C
	
	
return 0;
}
