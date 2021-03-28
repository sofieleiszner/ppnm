#include<stdio.h>
#include<math.h>
double myexp(double);
int main(){
	for(double x=1./8;x<=10;x+=1./8)
		printf("%g %g %g\n",x,myexp(x),exp(x));
return 0;
}



