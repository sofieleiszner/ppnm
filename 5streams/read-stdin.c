#include <stdio.h>
#include <math.h>

int main(){
	double x;
	int items;
	do {
	items = fscanf(stdin, "%lg", &x);
	printf("x = %g sin(x) = %g cos(x) = %g \n", x, sin(x), cos(x));
	} while(items != EOF);

return 0; 
}
