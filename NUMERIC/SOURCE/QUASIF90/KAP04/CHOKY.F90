![  {The Cholesky Decomposition}                                        
![  {The Cholesky Decomposition}*)                                      
      SUBROUTINE CHOKY (N, A, LDA, RS, X, Z, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Solving a linear system of equations                       *      
!                    A * X = RS                                  *      
!     for a symmetric, positive definite matrix A the            *      
!     Cholesky-method.                                           *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N    : order of the square matrix A                        *      
!     A    : 2-dimensional array A(1:LDA,1:N) containing the     *      
!            symmetric, positive definite matrix A. It is        *      
!            sufficient to specify only the elements of the      *      
!            upper triangle of A (including the elements of the  *      
!            main diagonal); only the elemente in A's upper      *      
!            triangle will be processed.                         *      
!     LDA  : leading dimension of A as defined in the calling    *      
!            program                                             *      
!     RS   : N-vector RS(1:N) containing the right hand side     *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     A    : 2-dimensional array A(1:LDA,1:N) which contains     *      
!            the Cholesky factors of A = R(TRANSP) * D * R       *      
!            for a diagonal D and a unit upper triangular matrix *      
!            R. The elements of R, excluding the diagonal ones,  *      
!            are stored in the upper triangle of A. The elements *      
!            of D appear on the main diagonal of A.              *      
!     X    : N-vector X(1:N) containing the solution of the      *      
!            system of equations                                 *      
!     MARK : error parameter                                     *      
!            MARK= 1 : ok                                        *      
!            MARK= 0 : numerically the matrix A is not stongly   *      
!                      nonsingular                               *      
!            MARK=-1 : A is not positive definite                *      
!                                                                *      
!     NOTE : If MARK = 1 , then the determinant of A can be      *      
!            calculated as follows:                              *      
!               DET A = A(1,1) * A(2,2) * ... * A(N,N)           *      
!                                                                *      
!                                                                *      
!     AUXILIARY PARAMETER:                                       *      
!     ====================                                       *      
!     Z    : N-vector Z(1:N)                                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: CHOKYP, CHOKYS, MACHPD                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author    : Gisela Engeln-Muellges                            *      
!  date      : 25.04.1988                                        *      
!  source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION A (LDA, N), RS (N), X (N), Z (N) 
!                                                                       
!  Cholesky factorization of the matrix A                               
!                                                                       
      CALL CHOKYP (N, A, LDA, Z, MARK) 
!                                                                       
!  Updating and backsustitution                                         
!                                                                       
      IF (MARK.EQ.1) THEN 
         CALL CHOKYS (N, A, LDA, RS, X, Z) 
      ENDIF 
      RETURN 
      END SUBROUTINE CHOKY                          
!                                                                       
!                                                                       
      SUBROUTINE CHOKYP (N, A, LDA, Z, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Factoring a symmetric, positive definite matrix A into     *      
!             A = R(TRANSP) * D * R  using the  Cholesky method. *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N    : order of the matrix A                               *      
!     A    : 2-dimensional array A(1:LDA,1:N) containing the     *      
!            symmetric, positive definite matrix A. It is        *      
!            sufficient to store only the elements of the upper  *      
!            triangle of A including the main diagonal due to    *      
!            symmetry.                                           *      
!     LDA  : leading dimension of A as defined in the calling    *      
!            program                                             *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     A    : 2-dimensional array A(1:LDA,1:N) containing the     *      
!            Cholesky Factors R and D in the array A where       *      
!            A = R(TRANSP) * D * R. The elements of R except for *      
!            the unit diagonal are stored in the upper triangle  *      
!            of A. The elements of D are on the main diagonal of *      
!            A.                                                  *      
!     MARK : error parameter                                     *      
!            MARK= 1 : ok                                        *      
!            MARK= 0 : numerically the matrix A is not strongly  *      
!                      nonsingular                               *      
!            MARK=-1 : A is not positive definite                *      
!                                                                *      
!                                                                *      
!     AUXILIARY PARAMETER:                                       *      
!     ====================                                       *      
!     Z    : N-vector Z(1:N)                                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author    : Gisela Engeln-Muellges                            *      
!  date      : 25.04.1988                                        *      
!  source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION A (LDA, N), Z (N) 
      MARK = 1 
!                                                                       
!  calculating the machine constant                                     
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = FMACHP * 2.0D0 
!                                                                       
!  determining the relative error bound                                 
!                                                                       
      EPS = 4.0D0 * FMACHP 
!                                                                       
!  calculating the absolute row sums of the matrix A,                   
!  with storage in the auxiliary vector Z                               
!                                                                       
      DO 20 I = 1, N 
         S = 0.0D0 
         DO 30 K = I, N 
            S = S + DABS (A (I, K) ) 
   30    END DO 
         DO 40 K = 1, I - 1 
            S = S + DABS (A (K, I) ) 
   40    END DO 
         IF (S.EQ.0.0D0) THEN 
            MARK = 0 
            RETURN 
         ENDIF 
         Z (I) = 1.0D0 / S 
   20 END DO 
!                                                                       
!  Factoring the matrix A                                               
!                                                                       
      DO 50 J = 1, N 
         DO 60 I = 1, J - 1 
            H = A (I, J) 
            A (I, J) = H / A (I, I) 
            DO 60 K = I + 1, J 
               A (K, J) = A (K, J) - H * A (I, K) 
   60    CONTINUE 
         IF (A (J, J) .LT.0.0D0) THEN 
            MARK = - 1 
            RETURN 
         ELSEIF (DABS (A (J, J) ) * Z (J) .LE.EPS) THEN 
            MARK = 0 
            RETURN 
         ENDIF 
   50 END DO 
      RETURN 
      END SUBROUTINE CHOKYP                         
!                                                                       
!                                                                       
      SUBROUTINE CHOKYS (N, A, LDA, RS, X, Z) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Finding the solution of the linear system of equations     *      
!                           A * X = RS                           *      
!     for a symmetric, positive definite matrix A that is        *      
!     given by its Cholesky factors by SUBROUTINE CHOKYP as      *      
!                   A = R(TRANSP) * D * R.                       *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N    : order of the matrix A                               *      
!     A    : output array of SUBROUTINE CHOKYP                   *      
!     LDA  : leading dimension of A as defined in the calling    *      
!            program                                             *      
!     RS   : N-vector RS(1:N) containing the right hand side     *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     X    : N-vector X(1:N) containing the solution of the      *      
!            system of equations                                 *      
!     RS   : overwritten with intermediate values                *      
!                                                                *      
!     N, A, LDA are not changed by this program                  *      
!                                                                *      
!                                                                *      
!     AUXILIARY PARAMETER:                                       *      
!     ====================                                       *      
!     Z    : N-vector Z(1:N)                                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author    : Gisela Engeln-Muellges                            *      
!  date      : 25.04.1988                                        *      
!  source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION A (1:LDA, 1:N), RS (1:N), X (1:N), Z (1:N) 
!                                                                       
!  updating                                                             
!                                                                       
      Z (1) = RS (1) 
      RS (1) = Z (1) / A (1, 1) 
      DO 10 J = 2, N, 1 
         Z (J) = RS (J) 
         DO 20 I = 1, J - 1, 1 
            Z (J) = Z (J) - A (I, J) * Z (I) 
   20    END DO 
         RS (J) = Z (J) / A (J, J) 
   10 END DO 
!                                                                       
!  backsubstitution                                                     
!                                                                       
      X (N) = RS (N) 
      DO 30 I = N - 1, 1, - 1 
         X (I) = RS (I) 
         DO 30 J = I + 1, N, 1 
            X (I) = X (I) - A (I, J) * X (J) 
   30 CONTINUE 
      RETURN 
      END SUBROUTINE CHOKYS                         
