#include<math.h>
#include<gsl/gsl_vector.h>
#include<assert.h>

int binsearch(int n,gsl_vector *x, double z){/* locates the interval for z by bisection */ 
	assert(gsl_vector_get(x,0)<=z && z<=gsl_vector_get(x,n-1));
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>gsl_vector_get(x,mid)) i=mid; else j=mid;
		}
return i;
}
