#include<gsl/gsl_matrix.h>
#include <stdio.h>
#include<gsl/gsl_blas.h>


//Functions from functions.c
void GenerateMatrix(gsl_matrix* A, int M, int N); 
void GKLBidiag(gsl_matrix* A, gsl_matrix* U, gsl_matrix* V, gsl_matrix* B); 
void matrix_print(char s[], gsl_matrix* M);
double detfromB(gsl_matrix* B); 
double detfromA(gsl_matrix* A); 
void SolveGKL(gsl_matrix* U, gsl_matrix*  B, gsl_matrix* Vt, gsl_vector* b, gsl_vector* y); 
void vector_print(char s[], gsl_vector* v);
void GenerateVector(gsl_vector* v); 
void InversGKL(gsl_matrix* U, gsl_matrix* B, gsl_matrix* V, gsl_matrix* Ainv);

int main(){
    printf("---------------------------------------------------------------------------------------------\n"); 

    printf("Exam, task 19: Golub-Kahan-Lanczos bidiagonalization\n"); 
    printf("Sofie Stampe Leiszner, 201705219\n"); 
    printf("---------------------------------------------------------------------------------------------\n"); 

    
    int N = 5; 
    int M = 5; 
    
    
    //Allocates memory
    gsl_matrix* A = gsl_matrix_alloc(M,N); 
    gsl_matrix* Ainv = gsl_matrix_alloc(N,N); 
    gsl_matrix* AAinv = gsl_matrix_alloc(N,N); 
    gsl_matrix* V = gsl_matrix_alloc(N,N); 
    gsl_matrix* U = gsl_matrix_alloc(M,M); 
    gsl_matrix* B = gsl_matrix_alloc(M,N); 
    gsl_matrix* BVt = gsl_matrix_alloc(M,N); 
    gsl_matrix* UBVt = gsl_matrix_alloc(M,N);  
    gsl_matrix* UtU = gsl_matrix_alloc(M,M); 
    gsl_matrix* UUt = gsl_matrix_alloc(M,M); 
    gsl_matrix* VtV = gsl_matrix_alloc(N,N); 
    gsl_matrix* VVt = gsl_matrix_alloc(N,N); 
    
    
    printf("Generates a random matrix with dimensions %d x %d\n", M, N); 
    GenerateMatrix(A, M, N); 
    matrix_print("A=",A); 
    
    printf("Applies the Golub-Kahan-Lanczos bidiagonalization algoritm\n"); 
    printf("The algoritm is implemented in the function GKLBidiag in functions.c\n"); 
    GKLBidiag(A, U, V, B); 
    printf("U, V and B are calculated\n");
    matrix_print("U = ", U); 
    matrix_print("V = ", V); 
    matrix_print("B = ", B); 
    printf("From the print of B, it is seen that it is upper bidiagonal\n"); 

        
    printf("Checks that U and V are orthogonal: Requires that U^T*U = U*U^T = I and V^T*V = V*V^T = I\n"); 
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, U, U, 0.0, UtU); //U^T*U
    gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, U, U, 0.0, UUt); //U*U^T
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, V, V, 0.0, VtV); //V^T*V
    gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, V, V, 0.0, VVt); //V*V^T
    matrix_print("U^T*U = ", UtU);
    matrix_print("U*U^T = ", UUt);
    matrix_print("V^T*V = ", VtV);
    matrix_print("V*V^T = ", VVt);


    printf("Check: that U*B*V^T = A\n"); 
    gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, B, V, 0.0, BVt); //BV^T
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, U, BVt, 0.0, UBVt); //U*BV^T
    matrix_print("Original A =", A); 
    matrix_print("U*B*V^T (should be A) =", UBVt); 
    
    printf("---------------------------------------------------------------------------------------------\n"); 
    //Linear equations
    printf("The algoritm can also be used to solve linear equations, Ax = b.\n"); 
    printf("Since A = U*B*V^T, then the the system of linear equations can be written as U*B*V^T*x=b .\n"); 
    printf("Since U is ortogonal: B*V^T*x=U^T*b\n"); 
    printf("Defines: y = V^T*x and solves B*y=U^T*b.\n"); 
    printf("And then calculates x = V*y.\n"); 
    printf("The function to do this is SolveGKL, which is implementes in functions.c\n"); 
    
    gsl_vector* b = gsl_vector_alloc(M); 
    gsl_vector* x = gsl_vector_alloc(M); 
    gsl_vector* Ax = gsl_vector_alloc(M); 
    printf("A random vector, b, is generated:"); 
    GenerateVector(b);
    vector_print("b = ", b); 
    printf("A*x=b is solved for x\n");
    SolveGKL(U, B, V, b, x);  
    vector_print("x = ", x);     
    printf("Checks that A * x = b\n"); 
    gsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 0.0, Ax); 
    vector_print("A*x =", Ax); 
    
    //Determinant
    printf("---------------------------------------------------------------------------------------------\n"); 
    printf("The algoritm can also be used to find the determinant of the matrix A\n"); 
    printf("Since A = U*B*V^T, then det(A)=det(U)det(B)det(V^T).\n"); 
    printf("Since U and V are orthogonale: U^2=I and V^2=I\n"); 
    printf("Therefore |det(A)|=|det(B)|\n"); 
    printf("Since B is triangular, det(B) = Πi Bii |\n");  

    printf("Calculates |det(A)|...\n"); 
    printf("|det(A)| = %g\n", detfromA(A)); 

    //Inverse
    printf("---------------------------------------------------------------------------------------------\n"); 
    printf("The algoritm can also be used to find the inverse of the matrix A\n"); 
    printf("The method is very similar to the one in Homework 2\n"); 
    printf("Simply solves A*x_k =e_k for k = 1,2,...,n with the linear equation solver.\n"); 
    printf("Where e_k are unit vectors with a 1 on the k'th element and zero everywhere else.\n"); 
    
    printf("Calculating the inverse of A...\n"); 
    InversGKL(U, B, V, Ainv);
    matrix_print("A^-1 = ", Ainv);
    printf("Checks that A*A^-1 = A^-1*A = I:\n");
    gsl_blas_dgemm(CblasTrans,CblasTrans, 1.0, A, Ainv, 0.0, AAinv); //V^T*V
    matrix_print("A*A^-1 =", AAinv); 
    gsl_blas_dgemm(CblasTrans,CblasTrans, 1.0, Ainv, A, 0.0, AAinv); //V^T*V
    matrix_print("A^-1*A =", AAinv); 
    printf("---------------------------------------------------------------------------------------------\n"); 
    
    gsl_matrix_free(A    ); 
    gsl_matrix_free(Ainv ); 
    gsl_matrix_free(AAinv); 
    gsl_matrix_free(V    ); 
    gsl_matrix_free(U    ); 
    gsl_matrix_free(B    ); 
    gsl_matrix_free(BVt  ); 
    gsl_matrix_free(UBVt );  
    gsl_matrix_free(UtU  ); 
    gsl_matrix_free(UUt  ); 
    gsl_matrix_free(VtV  ); 
    gsl_matrix_free(VVt  ); 
    
    return 0; 
}