![KA{P 8}{Linear and Nonlinear Approximation}                           
![  {Linear and Nonlinear Approximation}*)                              
![  {Normal Equations for Discrete Linear Least Squares}                
![  {Normal Equations for Discrete Linear Least Squares}*)              
      SUBROUTINE GADESM (M, X, F, W, LDA, N, C, A, B, Y, Z, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  SUBROUTINE GADESM determines the coefficients of a polynomial *      
!  of degree N that approximates a function f at the given nodes *      
!  in the discrete Gaussian least squares sense.                 *      
!  The linear system of normal equations is solved using the     *      
!  Cholesky method.                                              *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  M        : index of the last node.                            *      
!  X        : (N+1)-vector X(0:M) containing the nodes.          *      
!  F        : (N+1)-vector F(0:M) containing the function values *      
!             at the nodes.                                      *      
!  W        : (N+1)-vector W(0:M) containing the weights.        *      
!  LDA      : leading dimension of auxiliary matrix A as defined *      
!             in the calling program, LDA >= N+1.                *      
!  N        : degree of the approximating polynomial,            *      
!             2 <= N <= M.                                       *      
!                                                                *      
!                                                                *      
!  OUPUT PARAMETERS:                                             *      
!  =================                                             *      
!  C        : (N+1)-vector C(0:N) containing coefficients of the *      
!             approximating polynomial.                          *      
!  IERR     : = 0, no error.                                     *      
!             = 1, incorrect input parameter.                    *      
!             = 2, error in SUBROUTINE CHOKY.                    *      
!                                                                *      
!                                                                *      
!  HELP PARAMETER:                                               *      
!  ===============                                               *      
!  A        : 2-dim. array A(1:LDA,1:N+1).                       *      
!  B        : (N+1)-vector B(1:N+1).                             *      
!  Y        : (N+1)-vector Y(1:N+1).                             *      
!  Z        : (N+1)-vector Z(1:N+1).                             *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: CHOKY                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Guido Dubois                                       *      
!  date     : 05.30.87                                           *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!  declarations.                                                        
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (0:M), F (0:M), W (0:M), C (0:N), A (1:LDA, 1:N + 1), &
      B (1:N + 1), Y (1:N + 1), Z (1:N + 1)                             
!                                                                       
!  testing the input parameter.                                         
!                                                                       
      IERR = 1 
      IF (N.LT.2.OR.N + 1.GT.LDA.OR.M.LT.N) RETURN 
      IERR = 0 
!                                                                       
!  compute the first column of the system matrix for the normal equation
!  and the right-hand side.                                             
!                                                                       
      A (1, 1) = 0.0D0 
      B (1) = 0.0D0 
      DO 10 I = 0, M 
         A (1, 1) = A (1, 1) + W (I) 
         B (1) = B (1) + W (I) * F (I) 
   10 END DO 
      DO 20 J = 1, N 
         J1 = J + 1 
         A (J1, 1) = 0.0D0 
         B (J1) = 0.0D0 
         DO 30 I = 0, M 
            DUMMY = W (I) * X (I) **J 
            A (J1, 1) = A (J1, 1) + DUMMY 
            B (J1) = B (J1) + DUMMY * F (I) 
   30    END DO 
   20 END DO 
!                                                                       
!  compute the last row.                                                
!                                                                       
      DO 40 K = 1, N 
         K1 = K + 1 
         L = K + N 
         A (J1, K1) = 0.0D0 
         DO 50 I = 0, M 
            A (J1, K1) = A (J1, K1) + W (I) * X (I) **L 
   50    END DO 
   40 END DO 
!                                                                       
!  complete the matrix.                                                 
!                                                                       
      DO 60 K = 1, N 
         DO 70 I = 1, N 
            A (I, K + 1) = A (I + 1, K) 
   70    END DO 
   60 END DO 
!                                                                       
!  solve the system of normal equations. (The system matrix is          
!  positive definite, or MARK = 1 after CHOKY).                         
!                                                                       
      CALL CHOKY (N + 1, A, LDA, B, Y, Z, MARK) 
      IF (MARK.EQ.1) THEN 
         DO 80 J = 0, N 
            C (J) = Y (J + 1) 
   80    END DO 
      ELSE 
         IERR = 2 
      ENDIF 
      RETURN 
      END SUBROUTINE GADESM                         
