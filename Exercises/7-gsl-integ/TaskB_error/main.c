#include<stdio.h>
#include<math.h>
#include <gsl/gsl_integration.h>

double myerror(double z);


double f (double t, void * params) {
  double f = 2/sqrt(M_PI)*exp(-pow(t,2));
  return f;
}



int main(){
    printf("Task B\n");
    FILE * file = fopen("data.txt", "w");
	for(double x=-4;x<=4;x+=0.1/8){
		fprintf(file,"%10g %10g\n",x,myerror(x));}
    fclose(file);
    printf("Data is saved in data.txt.\n");
    printf("See erf.pyxplot.png for plot.\n");
return 0;
}

