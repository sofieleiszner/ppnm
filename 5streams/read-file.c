#include <math.h>
#include <stdio.h>

int main(int argc, char** argv) {
	FILE* my_in_stream = fopen(argv[1], "r"); 
	FILE* my_out_stream = fopen(argv[2], "w");
	int items;
	double x; 
	do {
	items = fscanf(my_in_stream, "%lg", &x); 
	fprintf(my_out_stream, "x = %g, cos(x) = %g \n", x, cos(x));
	} while(items != EOF);
	fclose(my_in_stream);
	fclose(my_out_stream);
} 
