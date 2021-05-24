#include <stdio.h>
#include <math.h>
#include<stdlib.h>



int main(){
    printf("Reads input from standard input\n");
	double x;
	int items;
	do {
	items = fscanf(stdin, "%lg", &x);
	fprintf(stdout, "x = %g sin(x) = %g cos(x) = %g \n", x, sin(x), cos(x));
	
    } while(items != EOF);

return 0; 
}
