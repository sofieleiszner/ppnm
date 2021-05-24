#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>

double Gamma(double);


int main(){
printf("Task A3 \n");
double xmax = 5, xmin = 1;
FILE * file = fopen("data.txt", "w"); 

for (double x = xmin; x<=xmax; x+=1.0/8) {
    fprintf(file, "%10g %10g %10g %10g\n" ,x, tgamma(x), gsl_sf_gamma(x), Gamma(x));
	}
fclose(file);
printf("Data is in data.txt \n");
printf("See plot in gamma.pyxplot.png \n");
return 0; 
}
