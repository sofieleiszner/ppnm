#include <stdio.h>
#include <limits.h>
#include <float.h>
#include "equal.c"


// Man kan ramme limit ved at compile med clang og køre uden -O option. 


int main(){
	//int i=1; while(i+1>i) {i++; if(i>100000) break;}
	//printf("my max int = %i \n", i);
	//i=1; 
	//for (i = 1; i+1>i; i++) {if(i>100000) break;};
	//printf("my max int = %i\n", i);
	//i=1; 	
	//do{
	//i++;
	//if(i>100000) break;} 
	//while(i+1>i);
	//printf("my max int = %i\n", i);
	//double x = 1; while(1+x != 1) {x/=2;} x*=2;
	//printf("x = %g\n", x);
	//float y  = 1; while(1+y != 1) {y/=2;} y*=2; 
	//printf("y = %f\n", y); 
	//long double z = 1; while(1+z != 1) {z/=2;} z*=2;
	//printf("z = %Lg\n", z);  
	int max = 100; // 2147483647/2;  skal rigtig være dette, men for speed 
	float sum_up_float = 0; 
	int i = 1;
	for (i = 1; i<max; i++) {
		sum_up_float += 1.0f/i;
	}
	printf("i = %d \n",i); 
	printf("sum_up_float = %g \n", sum_up_float);
	float sum_down_float = 1.0f/max;
	for (i = 1; i<max; i++) {
	sum_down_float +=1.0f/(max-i);
	}
	printf("i = %d \n", i);
	printf("sum_down_float = %g \n", sum_down_float); 
	// Sum up float lægges de store tal sammen først, det giver: 15.4
	// Sum down float lægges de små tal sammen først, det giver 18.8
	// Det giver højere præcision at lægge de små tal sammen først
	//For hvis man først har lagt mange store tal sammen og derefter lægge 
	//de små tal til, så bliver de bare decimaler, som ikke kan tages med
	// pga. bit begrænsning. 
	double sum_up_double = 0; 
	for (i = 1; i<max; i++) {
		sum_up_double += 1.0f/i;
	}
	printf("i = %d \n",i); 
	printf("sum_up_double = %g \n", sum_up_double);
	double sum_down_double = 1.0f/max;
	for (i = 1; i<max; i++) {
	sum_down_double +=1.0f/(max-i);
	}
	printf("i = %d \n", i);
	printf("sum_down_double = %g \n", sum_down_double); 
	int equal_result = equal(8.3, 8.1, 0.5, 0.1);
	printf("equal_result = %d \n", equal_result);
return 0;
}
