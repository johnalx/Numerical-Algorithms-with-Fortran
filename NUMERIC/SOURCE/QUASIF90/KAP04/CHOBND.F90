![          {Linear Systems with Band Matrices using Pivots}*)          
      SUBROUTINE CHOBND (N, M, AP, RS, X, IFLAG, Z) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  CHOBND solves a linear system of equations                    *      
!                 A * X = RS                                     *      
!  for a symmetric, positive definite banded matrix A in         *      
!  condensed form using the Cholesky method.                     *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N      : number of rows of A                                  *      
!  M      : number of nonzero codiagonals of A                   *      
!  AP     : array AP(1:N,1:M+1), containing the upper triangular *      
!           entries of A in condensed form.                      *      
!           The condensation of A is achieved by shifting the    *      
!           rows of A successively to the left until the diagonal*      
!           and upper codiagonals of A appear as columns.        *      
!           The main diagonal of A forms the first column of the *      
!           condensed matrix; the first codiagonal is in its     *      
!           second column, its last element is set equal to zero;*      
!           etc. until the (M+1)-st column of the condensed      *      
!           matrix contains the M-th codiagonal of A in positions*      
!           1 to N-M and zeros below.                            *      
!           The following code will condense a matrix A into AP: *      
!                  DO 10 I=1,N                                   *      
!                     DO 20 K=I,MIN(N,I+M)                       *      
!                        AP(I,K-I+1) = A(I,K)                    *      
!               20    CONTINUE                                   *      
!               10 CONTINUE                                      *      
!                                                                *      
!  RS     : N-vector RS(1:N), the right hand side of the system  *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  AP     : Same as A, but overwritten with the Cholesky factors *      
!           R(TRANS)*D*R in condensed form: R is an upper band   *      
!           triangular matrix with unit diagonal and M bands,    *      
!           D is a diagonal matrix.                              *      
!           The diagonal entries of D form the first column of   *      
!           AP, successive columns of AP contain the codiagonals *      
!           of R. The unit diagonal of R is not stored and       *      
!           neither is R(TRANS).                                 *      
!  X      : N-vector X(1:N), the solution vector                 *      
!  IFLAG  :  error parameter:                                    *      
!            1: all is ok                                        *      
!            0: Matrix is numerically singular                   *      
!           -1: Matrix is numerically not positive definite      *      
!                                                                *      
!           If  IFLAG = 1, the  determinant of A is given by:    *      
!              DET(A) = AP(1,1) * AP(2,1) * ... * AP(N,1)        *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETER:                                          *      
!  ====================                                          *      
!  Z      : N-vector Z(1:N)                                      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines:                                         *      
!                                                                *      
!  CHOBDZ  Cholesky factorization of AP                          *      
!  CHOBDL  Updating and backsustitution                          *      
!  MACHPD  machine constant                                      *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author     : Elmar Pohl                                       *      
!  Date       : 11.15.1991                                       *      
!  Sourcee    : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      DIMENSION AP (1:N, 1:M + 1), RS (1:N), X (1:N), Z (1:N) 
!                                                                       
!     Factor  A = R(TRANS)*D*R                                          
!                                                                       
      CALL CHOBDZ (N, M, AP, IFLAG, Z) 
!                                                                       
!     Updating and backsubstitution                                     
!                                                                       
      IF (IFLAG.EQ.1) CALL CHOBDL (N, M, AP, RS, X, Z) 
!                                                                       
      RETURN 
      END SUBROUTINE CHOBND                         
!                                                                       
!                                                                       
      SUBROUTINE CHOBDL (N, M, AP, RS, X, Z) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  CHOBDL solves a linear system of equations                    *      
!                 A * X = RS                                     *      
!  for a symmetric, positive definite banded matrix A in         *      
!  condensed form using the Cholesky method.                     *      
!  Here the matrix A has been overwritten in AP by the subroutine*      
!  CHOBDZ with its Cholesky factors  A = R(TRANS)*D*R  in        *      
!  condensed form.                                               *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N      : number of rows of AP                                 *      
!  M      : number of nonzero codiagonals of A                   *      
!  AP     : array AP(1:N,1:M+1), the output of SUBROUTINE CHOBDZ *      
!           with the Cholesky factors of A in condensed form.    *      
!  RS     : N-vector RS(1:N), the right hand side of the system  *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!                                                                *      
!  X      : N-vector X(1:N), the solution vector                 *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETER:                                          *      
!  ====================                                          *      
!  Z      : N-vector Z(1:N)                                      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author     : Elmar Pohl                                       *      
!  Date       : 11.15.1991                                       *      
!  Sourcee    : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      DIMENSION AP (1:N, 1:M + 1), RS (1:N), X (1:N), Z (1:N) 
!                                                                       
!     Solve  R(TRANS) * Z = RS                                          
!                                                                       
      DO 10 J = 1, N 
         Z (J) = RS (J) 
         DO 20 I = MAX0 (1, J - M), J - 1 
            Z (J) = Z (J) - AP (I, J - I + 1) * Z (I) 
   20    END DO 
   10 END DO 
!                                                                       
!     Solve  R * X = D^-1 * Z                                           
!                                                                       
      DO 30 I = N, 1, - 1 
         X (I) = Z (I) / AP (I, 1) 
         DO 40 J = I + 1, MIN0 (N, I + M) 
            X (I) = X (I) - AP (I, J - I + 1) * X (J) 
   40    END DO 
   30 END DO 
      RETURN 
      END SUBROUTINE CHOBDL                         
!                                                                       
!                                                                       
      SUBROUTINE CHOBDZ (N, M, AP, IFLAG, Z) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  CHOBDZ factors a symmetric, positive definite banded matrix AP*      
!  given in condensed form into  R(TRANS)*D*R  using the         *      
!  Cholesky method. Here D is a diagonal matrix and R is a unit  *      
!  diagonal upper tringular banded matrix with as many           *      
!  codiagonals as the original A. The output is again stored in  *      
!  condensed form with D in the first column and the codiagonal  *      
!  bands of R in subsequent columns.                             *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N      : number of rows of AP                                 *      
!  M      : number of codiagonals in A                           *      
!  AP     : array AP(1:N,1:M+1), containing the upper triangle   *      
!           of A in condensed form. For a code that achieves this*      
!           for a given symmetric matrix, see CHOBND.            *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  AP     : Same as A, but overwritten with the Cholesky factors *      
!           R(TRANS)*D*R in condensed form: R is an upper band   *      
!           triangular matrix with unit diagonal and M bands,    *      
!           D is a diagonal matrix.                              *      
!           The diagonal entries of D form the first column of   *      
!           AP, successive columns of AP contain the codiagonals *      
!           of R. The unit diagonal of R is not stored and       *      
!           neither is R(TRANS).                                 *      
!  X      : N-vector X(1:N), the solution vector                 *      
!  IFLAG  :  error parameter:                                    *      
!            1: all is ok                                        *      
!            0: Matrix is numerically singular                   *      
!           -1: Matrix is numerically not positive definite      *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETERS:                                         *      
!  ====================                                          *      
!  Z      : N-vector Z(1:N)                                      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines:                                         *      
!                                                                *      
!  MACHPD  Machine constant                                      *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Elmar Pohl                                        *      
!  Date      : 11.15.1991                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      DIMENSION AP (1:N, 1:M + 1), Z (1:N) 
!                                                                       
!     Find the machine constant                                         
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = FMACHP * 2.0D0 
!                                                                       
!     Set a bound for the relative accuracy                             
!                                                                       
      EPS = 4.0D0 * FMACHP 
!                                                                       
!     Store row sum moduli of each row in Z,                            
!     use Z to check nonsingularity later                               
!                                                                       
      DO 20 I = 1, N 
         S = 0.0D0 
         DO 30 K = I, MIN0 (N, I + M) 
            S = S + DABS (AP (I, K - I + 1) ) 
   30    END DO 
         DO 40 K = MAX0 (1, I - M), I - 1 
            S = S + DABS (AP (K, I - K + 1) ) 
   40    END DO 
         IF (S.EQ.0.0D0) THEN 
            IFLAG = 0 
            RETURN 
         ENDIF 
         Z (I) = S 
   20 END DO 
!                                                                       
!     Cholesky decomposition                                            
!                                                                       
      DO 50 J = 1, N 
         DO 60 I = MAX0 (1, J - M), J - 1 
            H = AP (I, J - I + 1) 
            AP (I, J - I + 1) = H / AP (I, 1) 
            DO 70 K = I + 1, J 
               AP (K, J - K + 1) = AP (K, J - K + 1) - H * AP (I, K - I &
               + 1)                                                     
   70       END DO 
   60    END DO 
         IF (AP (J, 1) .LT.0.0D0) THEN 
            IFLAG = - 1 
            RETURN 
         ELSEIF (DABS (AP (J, 1) ) / Z (J) .LE.EPS) THEN 
            IFLAG = 0 
            RETURN 
         ENDIF 
   50 END DO 
!                                                                       
      IFLAG = 1 
      RETURN 
      END SUBROUTINE CHOBDZ                         
