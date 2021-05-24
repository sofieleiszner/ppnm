#include <math.h>
#include <stdio.h>

int main(int argc, char** argv) {

	FILE* my_in_stream = fopen(argv[1], "r"); 
	FILE* my_out_stream = fopen(argv[2], "w");
	int items;
	double x; 
	while((items = fscanf(my_in_stream, "%lg", &x)) != EOF) {
	fprintf(my_out_stream, "x = %g, cos(x) = %g \n", x, cos(x));
	};
	fclose(my_in_stream);
	fclose(my_out_stream);
    return 0; 
} 
