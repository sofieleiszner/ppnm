#include<stdio.h>
#include<math.h>


double f (double x, void * params) {
  double f = log(x)/sqrt(x);
  return f;
}
double myQAGS();


int main(){
    printf("Task A\n");
    printf("Integral from 0 to 1 of ln(x)/sqrt(x) = %f\n",myQAGS());
return 0;
}