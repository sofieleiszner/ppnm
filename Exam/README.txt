Exam project: "Golub-Kahan-Lanczos bidiagonalization"
by Sofie Stampe Leiszner (201705219)

The implemented program can represent a real matrix, A, in the form A = UBV^T 
(where U and V are othrogonal matrices and B is a diagonal matrix). 

A random matrix is generated, and U, B and V are calculated from it. 
Several checks are made to ensure that the algorithm functions as intended. 

The implemented algorithm is then used to solve linear equations of the form Ax=b
By using A=UBV^T
And then rewriting Ax=b to UBV^Tx=b
And since U is ortogonal: BV^Tx=U^Tb
Introducing the variable y: y = V^Tx
And then solving the equation upper triangular system By=b for y by backsubstitution. 
Then since V is ortogonal: x = Vy

Additionally, the determinant is determined: 
det(A)=det(U)det(B)det(V^T)
Since U and V are ortogonal: |det(U)| = |det(V)| = 1
|det(A)|=|det(B)|
Since B is triangular: |det(A)| = |Πi Bii|

Finally, the algoritm is used to determinde the inverse of A
By solving A*x_k = e_k for k =1,2,...,n with the linear equation solver 
explained above

All the above methods are used on the randomly generated matrix, A. 
Several checks are made. 
See the out.txt file for more. 