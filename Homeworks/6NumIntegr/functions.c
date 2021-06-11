#include <math.h>
#include <stdio.h>
#include <assert.h>

double global_error;
//int calls;

double integrate(double fg(double f(double x), double x, double a), double f(double x), double a0, double a, double b, double acc, double eps, double f2, double f3, int nrec)
{
	assert(nrec<1e6);
	double f1=fg(f, a+(b-a)/6.0, a0), f4=fg(f, a+5.0*(b-a)/6.0, a0);
	double Q=(2.0*f1+f2+f3+2.0*f4)/6.0*(b-a), q=(f1+f4+f2+f3)/4.0*(b-a);
	double tolerance=acc+eps*fabs(Q), error=fabs(Q-q);
	if (error < tolerance) {
		global_error += pow(error,2);
		return Q;
	}
	else {return integrate(fg, f, a0, a,(a+b)/2.0,acc/sqrt(2.),eps, f1, f2, nrec+1)+
				integrate(fg, f,a0,(a+b)/2.0,b,acc/sqrt(2.),eps, f3, f4, nrec+1);
	}
}

double g0(double f(double x), double x, double a){
	return f(x);
}

double g1(double f(double x), double x, double a){
			return (f((1-x)/x) + f(-(1.0-x)/x))*1.0/pow(x,2);
		}

double g2(double f(double x), double x, double a){
			return f(a + (1.0 - x)/x)/(x*x);
		}

double g3(double f(double x), double x, double b){
			return f(b - (1.0-x)/x)*1.0/pow(x,2);
		}

double infintegrate(double f(double x), double a, double b, double acc, double eps){
	global_error = 0;
	if (isinf(a)==-1 && isinf(b)==1){
		double A = 0.0, B = 1.0;
		double f2=g1(f, A+2.0*(B-A)/6.0, 0), f3=g1(f, A+4.0*(B-A)/6.0, 0); int nrec=0;
		return integrate(g1,f, a, A, B, acc, eps, f2, f3, nrec);
	}
	else if (isinf(b)==1){
		double A = 0.0, B = 1.0;
		double f2=g2(f, A+2.0*(B-A)/6, a), f3=g2(f, A+4.0*(B-A)/6, a); int nrec=0;
		return integrate(g2,f, a, A, B, acc, eps, f2, f3, nrec);
	}
	else if (isinf(a)==-1){
		double A = 0.0, B = 1.0;
		double f2=g3(f, A+2.0*(B-A)/6.0, b), f3=g3(f, A+4.0*(B-A)/6.0, b); int nrec=0;
		return integrate(g3,f, b, A, B, acc, eps, f2, f3, nrec);
	}
	else {
		double f2=g0(f, a+2.0*(b-a)/6.0, 0), f3=g0(f, a+4.0*(b-a)/6.0, 0); int nrec=0;
		return integrate(g0,f, a, a, b, acc, eps, f2, f3, nrec);
	}
}