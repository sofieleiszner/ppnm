#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<assert.h>
#include <gsl/gsl_interp.h>

int binsearch(int n,gsl_vector *x,double Z);
typedef struct {gsl_vector *x, *y, *b, *c;} qspline;
typedef struct {gsl_vector *x, *y, *b, *c, *d;} cubic_spline;

qspline* qspline_alloc(gsl_vector *x, gsl_vector* y);
cubic_spline* cubic_spline_alloc(gsl_vector *x, gsl_vector* y);

double linterp(gsl_vector * x,gsl_vector * y, double Z);
double linterp_integ(gsl_vector *x,gsl_vector *y, double Z);

double qspline_eval(qspline *s, double z);
double quad_integ(qspline *s, double z);
double quad_diff(qspline *s,double z);


double cubic_spline_eval(cubic_spline*s, double z);
double cubic_integ(cubic_spline*s, double z);
double cubic_diff(cubic_spline*s, double z);

void quad_spline_free(qspline *s);
void cubic_spline_free(cubic_spline *s);

int main(){
//	int N=10;

//finder længden af filen hvori vores punkter ligger og skriver det i heltallet N.
	int N=0;
	FILE* fp=fopen("input.dat","r");
	int ch;
	while(!feof(fp))
		{
  		ch = fgetc(fp);
  		if(ch == '\n')
  			{
    			N++;
  			}
		}
	fclose(fp);

//Memory is allocated to two vector x and y, which will contain the x and y values for each point.
	gsl_vector * x = gsl_vector_alloc (N);
	gsl_vector * y = gsl_vector_alloc (N);


//loader punkter fra fil ind i vektorerner x og y.
	double i;
	double j;
	int items;
	int n=0;
	FILE* my_in_stream=fopen("input.dat","r");
	while((items=fscanf(my_in_stream,"%lg %lg",&i,&j))!=EOF){
		gsl_vector_set (x, n, i);
		gsl_vector_set (y, n, j);
		n++;
		}
	fclose(my_in_stream);

// GSL linear interpolator.
	double gslx[N],gsly[N];

	for(int i=0;i<N;i++){
		gslx[i]=gsl_vector_get(x,i);
		gsly[i]=gsl_vector_get(y,i);
	}
	gsl_interp* linear= gsl_interp_alloc(gsl_interp_linear,N);
	gsl_interp_init(linear,gslx,gsly,N);

	gsl_interp* cubic=gsl_interp_alloc(gsl_interp_cspline,N);
	gsl_interp_init(cubic,gslx,gsly,N);
    printf("Task A\n");
    printf("Data from the file input.dat is loaded into gsl_vectors x and y\n"); 
    printf("The linear interpolation function is in linterp.c\n"); 
    printf("Integration of the linear spline is implemented in linterp_integ.c\n"); 
    printf("See plot in fig-linterp-pyxplot.png\n");
    printf("Task B\n");
    printf("Quadratic spline is implemented in quad.c\n");
    printf("And the corresponding derivaive and integral in quad_diff and quad_integ\n");
    printf("See plot in fig-quadterp-pyxplot.png\n");
    printf("Task C\n");
    printf("Quadratic spline is implemented in cubic.c\n");
    printf("And the corresponding derivaive and integral in cubic_diff and cubic_integ\n");
    printf("See plot in fig-cubicterp-pyxplot.png\n");
    printf("All data for linear, quadratic and cubic spline is put into linterp.dat\n");
    printf("All implementations are compared with GSL spline routines\n");
	FILE * f=fopen("linterp.dat","w");
	double z=gsl_vector_get(x,0);
	qspline *s=qspline_alloc(x,y);
	cubic_spline *ss=cubic_spline_alloc(x,y);
	for((z=z);z<=N;z+=1./16){
		double yiz=linterp(x,y,z);
		double integraltilz=linterp_integ(x,y,z);
		double quadterp=qspline_eval(s,z);
		double quadterptilz=quad_integ(s,z);
		double quadterpdiffiz=quad_diff(s,z);
		double interp_l=gsl_interp_eval(linear,gslx,gsly,z,NULL);
		double integ_l=gsl_interp_eval_integ(linear,gslx,gsly,gslx[0],z,NULL);
		double cubicterp=cubic_spline_eval(ss,z);
		double cubicterpdiffiz=cubic_diff(ss,z);
		double cubicterptilz=cubic_integ(ss,z);
		double interp_cubic=gsl_interp_eval(cubic,gslx,gsly,z,NULL);
		double integ_cubic=gsl_interp_eval_integ(cubic,gslx,gsly,gslx[0],z,NULL);
		double diff_cubic=gsl_interp_eval_deriv(cubic,gslx,gsly,z,NULL);
		fprintf(f,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",z,yiz,integraltilz,quadterp,quadterptilz,quadterpdiffiz,interp_l,integ_l,cubicterp,cubicterpdiffiz,cubicterptilz,interp_cubic,integ_cubic,diff_cubic);
		}
	fclose(f);

//quad_spline_free(s);
//cubic_spline_free(ss);
gsl_interp_free(linear);
gsl_vector_free (x);
gsl_vector_free (y);

printf("Program done\n");

return 0;
}
