#include <stdio.h>
#include <math.h> 
#include <gsl/gsl_sf_erf.h>

double Erf(double); 

int main(){
    printf("Task A \n");
    printf("Task A1 \n");
    printf("Installation done \n");
    printf("Task A2 \n");
    printf("Error function approximation implemented in erf.c \n");
	double xmin = -2, xmax =2; 
    FILE * file = fopen("data.txt", "w"); 
	for(double x = xmin; x<=xmax; x+=1.0/8) {
	fprintf(file, "%10g %10g %10g %10g\n",x,erf(x),gsl_sf_erf(x),Erf(x));
	}
    fclose(file);
    printf("Data is in data.txt \n");
    printf("See plot in erf.pyxplot.png \n");
return 0;
}
