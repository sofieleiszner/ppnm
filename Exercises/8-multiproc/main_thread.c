#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <pthread.h>


typedef struct{
    double x; 
    double y;} 
    point; 
    
typedef struct{
    unsigned int seed; 
    double N;} 
    seedN; 

int inorout(point p){
    double x_vrdi = p.x;
    double y_vrdi = p.y;
    if(sqrt(pow(x_vrdi,2)+pow(y_vrdi,2))<1){return 1;}
    else {return 0;}
}

void* howmanyinside(void* arg) //N is number of points    
 {  
    seedN * param = (seedN*) arg; 
    int N = (*param).N; 
    unsigned int seed = (*param).seed;
    int sum = 0 ;
    for (int i = 0; i < N; i++) {
    point p1 = {.x=((double)rand_r(&seed)/RAND_MAX),.y=((double)rand_r(&seed)/RAND_MAX)};     
    sum += inorout(p1);
 }
    (*param).N = sum; 
return NULL;
}



int main() {
    printf("Task A\n");
    int N = 1000000;
    pthread_t threadx, thready;
    pthread_attr_t* attributes = NULL;
    seedN N1 = {.seed = 1, .N = N/2};
    seedN N2 = {.seed = 2, .N = N/2};
    pthread_create( &threadx, attributes, howmanyinside, (void*)&N1 );
	pthread_create( &thready, attributes, howmanyinside, (void*)&N2 );
	void* returnvalue = NULL;
	pthread_join(threadx,returnvalue);
	pthread_join(thready,returnvalue);
    double Nin = N1.N + N2.N;
    printf("Pi-approksimation = %g\n ", 4*((double)Nin/N));
return 0; 
}