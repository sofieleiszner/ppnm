
#include "komplex.h"
#include "stdio.h"
#define TINY 1e-6

int main(){
    printf("Task 1: See komplex.h and komplex.c \n"); 
    
    printf("Task 2: \n"); 
	komplex a = {1,2};
	komplex b = {3,4};
	printf("testing komplex_add...\n");
	komplex r = komplex_add(a,b);
	komplex R = {4,6};
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a+b should =", R);
	komplex_print("a+b actually =", r);
    
    printf("Task 3: See makefile \n"); 

return 0; 
}
