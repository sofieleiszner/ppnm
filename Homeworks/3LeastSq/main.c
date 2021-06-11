#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_blas.h>


void GenerateMatrix(gsl_matrix* A, int N, int M);
void matrix_print(char s[], gsl_matrix* M);
void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void vector_print(char s[], gsl_vector* v);
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);
void leastsq(gsl_vector* x, gsl_vector* y, gsl_vector* dy,gsl_vector* c,gsl_vector* dc, double f(int k, double z));
double f(int k, double z);

int main(){
    printf("TAsk A, B and C\n");
   // int N = 6; int M = 3;

    gsl_vector* t = gsl_vector_alloc(9);
    gsl_vector* y = gsl_vector_alloc(9);
    gsl_vector* dy = gsl_vector_alloc(9);
	gsl_vector* c = gsl_vector_alloc(2);
	gsl_vector* dc = gsl_vector_alloc(c->size);

	double i;
    double j;
	double k;
	double l;
	double m;
    int items;
    int n=0;
    printf("Reads data from data.txt \n");
    FILE* my_in_stream=fopen("data.txt","r");
    while((items=fscanf(my_in_stream,"%lg %lg %lg %lg %lg",&i,&j,&k,&l,&m))!=EOF){
        gsl_vector_set (t, n, i);
        gsl_vector_set (y, n, j);
		gsl_vector_set (dy, n, m);
        n++;
        }
    fclose(my_in_stream);

	gsl_vector* y_ln = gsl_vector_alloc(9);
	gsl_vector_memcpy(y_ln,y);
	for (int i=0; i<y_ln->size; i++){
		gsl_vector_set(y_ln,i,log(gsl_vector_get(y_ln,i)));
	}
	printf("Does a least square fit with the function leastsq\n");
	leastsq(t, y_ln, dy, c, dc, f);

    printf("Prints the fit into datafile.txt along with the two fits with the \n");
    printf("sum  in k over (ck+δck)*fk(x) and  sum in k over (ck-δck)*fk(x)\n");

    FILE * f1=fopen("datafile.txt","w");
	for (double x = 1.0/16; x<15; x+=1.0/8) {
		double sum = 0;
		double upsum = 0;
		double downsum = 0;
		for (int k = 0; k<c->size; k++) {
			sum += gsl_vector_get(c,k)*f(k,x);
			upsum += (gsl_vector_get(c,k)+gsl_vector_get(dc,k))*f(k,x);
			downsum += (gsl_vector_get(c,k)-gsl_vector_get(dc,k))*f(k,x);
		}

		fprintf(f1, "%g %g %g %g\n", x, sum, upsum, downsum);
	}


	printf("Half-life time = %g days\n", log(2)/gsl_vector_get(c,1));
	printf("Uncertainty in half-life time = %g days\n", log(2)/pow(gsl_vector_get(c,1),2)*gsl_vector_get(dc,1)); //fejlophobningslov
	printf("Half-life time, table value (today) = 3.6319 days\n"); //https://en.wikipedia.org/wiki/Isotopes_of_radium
	printf("Half-life time uncertainty (today) = 0.0023 days\n");
	printf("The half-life time does not agree with modern value within the uncertainty.\n");
	printf("But the values are not very different and the data is very old, so thats OK.\n");

	vector_print("Vector of coefficients, c = ", c);
	vector_print("Uncertanties of the coeficients, dc = ", dc);
    printf("Fit: ln(y)=(%.3g+-%.2g)+(%.3g+-%.1g)*t", gsl_vector_get(c,0), gsl_vector_get(dc,0), gsl_vector_get(c,1), gsl_vector_get(dc,1)); 
    printf("Plot of the three fits and datapoints with errorbars is shown in least.pyxplot.png\n"); 
	gsl_vector_free(t);
    gsl_vector_free(y);
    gsl_vector_free(dy);
    gsl_vector_free(c);
    gsl_vector_free(y_ln);

  //gsl_matrix* A = gsl_matrix_alloc(N, M);
  //gsl_matrix* A_save = gsl_matrix_alloc(N, M);
  //gsl_matrix* R = gsl_matrix_alloc(M, M);


  //  GenerateMatrix(A, N, M);
  //  gsl_matrix_memcpy(A_save,A);
  //
  //  matrix_print("Random matrix A = ", A);
  //
  //  GS_decomp(A, R);
  //
  //  printf("\nResults after decomposition: \n");
  //  matrix_print("R = ", R);
  //  matrix_print("Q = ", A);
  //
  //  printf("\nChecks:\n");
  //  //check that Q^T*Q;
  //  printf("Check that Q^T*Q = 1: \n");
  //  // gsl_blas_dgemm()
  //  gsl_matrix* QTQ = gsl_matrix_alloc(M, M);
  //  gsl_matrix* AT = gsl_matrix_alloc(M, N);
  //  gsl_matrix_transpose_memcpy(AT, A);
  //  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AT, A, 0.0, QTQ);
  //  matrix_print("Q^T*Q =", QTQ);
  //
  //  printf("Check that QR=A:\n");
  //  gsl_matrix* QR = gsl_matrix_alloc(N, M);
  //  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, R, 0.0, QR);
  //  matrix_print("QR=", QR);
  //
  //
  //  gsl_vector* b=gsl_vector_alloc(N);
	//for (int i=0;i<N;i++)
	//	gsl_vector_set(b,i,(double)rand()/RAND_MAX);
  //  gsl_vector* x = gsl_vector_alloc(M);
  //  GS_solve(A, R, b, x);
  //  //Checks Ax = b
  //  gsl_vector* Ax = gsl_vector_alloc(N);
  //  gsl_blas_dgemv(CblasNoTrans, 1.0, A_save, x, 0.0, Ax);
  //  matrix_print("A = ", A_save);
  //  vector_print("x=", x);
  //
  //  vector_print("Ax =", Ax);
  //  vector_print("b=", b);
  //
  //
  //
  //  gsl_matrix_free(A);
  //  gsl_matrix_free(R);
  //  gsl_matrix_free(QTQ);
  //  gsl_matrix_free(AT);
  //  gsl_matrix_free(QR);
  //  gsl_vector_free(b);
  //  gsl_vector_free(x);
  //  gsl_vector_free(Ax);
  //  gsl_matrix_free(A_save);


return 0;
}