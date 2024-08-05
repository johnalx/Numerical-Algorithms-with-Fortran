![KA{P 7}{Eigenvalues and Eigenvectors of Matrices}                     
![       {Eigenvalues and Eigenvectors of Matrices}*)                   
![  {Vector Iteration for the Dominant Eigenvalue}                      
![  {Vector Iteration for the Dominant Eigenvalue and the               
![   Associated Eigenvector of a Matrix}*)                              
      SUBROUTINE EIGVAL (A, N, LDA, M0, EPSI, X, Y, Z, EV, IERR) 
!                                                                       
!*****************************************************************      
!*                                                               *      
!*  Determining the eigenvalue of largest magnitude of an n by n *      
!*  matrix A with the corresponding eigenvector by vector        *      
!*  iteration                                                    *      
!*                                                               *      
!*                                                               *      
!*   INPUT PARAMETERS:                                           *      
!*   =================                                           *      
!*   A   : 2-dimensional DOUBLE PRECISION array A(1:LDA,1:N);    *      
!*         the input matrix                                      *      
!*   N   : order of A                                            *      
!*   LDA : leading dimension of A as defined in the calling      *      
!*         program                                               *      
!*   M0  : maximum iteration number                              *      
!*   EPSI: desired relative accuracy                             *      
!*         (larger than 1E-12) (DOUBLE PRECISION)                *      
!*   X   : N-vector X(1:n) in DOUBLE PRECISION                   *      
!*                                                               *      
!*                                                               *      
!*   OUTPUT PARAMETERS:                                          *      
!*   ==================                                          *      
!*   Y   : N-vector Y(1:n) in DOUBLE PRECISION; the eigenvector  *      
!*   Z   : N-vector Z(1:n) in DOUBLE PRECISION; the residual     *      
!*         vector A * Y - EV * Y                                 *      
!*   EV  : the dominant eigenvalue in DOUBLE PRECISION           *      
!*   IERR: error parameter:                                      *      
!*         =0: run was successfully completed                    *      
!*         =1: maximum number of iterations was reached, i.e.,   *      
!*             eigenvalue/vector is complex or the problem is    *      
!*             poorly conditioned                                *      
!*                                                               *      
!*---------------------------------------------------------------*      
!*                                                               *      
!*  subroutines required: DBNORM, MAVE, QUOT, QSCAL              *      
!*                                                               *      
!*                                                               *      
!*                                                               *      
!*****************************************************************      
!*                                                               *      
!*  author   : Juergen Beckmann                                  *      
!*  date     : 10.24.1985                                        *      
!*  source   : FORTRAN 77                                        *      
!*                                                               *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (LDA, LDA), Z (LDA), X (LDA), Y (LDA) 
      IERR = 0 
      DO 100 I = 1, N 
         Y (I) = 1.0D0 
  100 END DO 
      CALL DBNORM (Y, N) 
      EV = 0.0D0 
      ITER = 0 
  200 ITER = ITER + 1 
      EM = EV 
      DO 300 I = 1, N 
         X (I) = Y (I) 
  300 END DO 
      CALL MAVE (A, X, Y, N, LDA) 
      EV = QUOT (X, Y, N) 
      CALL DBNORM (Y, N) 
      IF (ITER.EQ.1) THEN 
         GOTO 200 
      ENDIF 
      DO 400 I = 1, N 
         Z (I) = Y (I) - X (I) 
  400 END DO 
      CALL QSCAL (Z, Z, S, N) 
      S = DSQRT (S) 
      IF (ITER.EQ.M0) THEN 
         IERR = 1 
         GOTO 500 
      ENDIF 
      IF (S.GT.EPSI.OR.DABS (EV - EM) .GT.EPSI) THEN 
         GOTO 200 
      ENDIF 
  500 RETURN 
      END SUBROUTINE EIGVAL                         
!                                                                       
!                                                                       
      SUBROUTINE MAVE (A, X, Y, N, LDA) 
!                                                                       
!*****************************************************************      
!*                                                               *      
!*  compute  Y = A * X                                           *      
!*                                                               *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (LDA, LDA), X (LDA), Y (LDA) 
      DO 200 I = 1, N 
         Y (I) = 0.0D0 
         DO 100 J = 1, N 
            Y (I) = Y (I) + A (I, J) * X (J) 
  100    END DO 
  200 END DO 
      RETURN 
      END SUBROUTINE MAVE                           
!                                                                       
!                                                                       
      SUBROUTINE DBNORM (X, N) 
!                                                                       
!*****************************************************************      
!*                                                               *      
!*  Normalizes a DOUBLE PRECISION vector to euclidean length 1   *      
!*                                                               *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (N) 
      CALL QSCAL (X, X, S, N) 
      S = DSQRT (S) 
      IF (S.NE.0.0D0) THEN 
         DO 100 I = 1, N 
            X (I) = X (I) / S 
  100    END DO 
      ENDIF 
      RETURN 
      END SUBROUTINE DBNORM                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION QUOT (X, Y, N) 
!                                                                       
!*****************************************************************      
!*                                                               *      
!*  auxiliary routine for EIGVAL                                 *      
!*                                                               *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (N), Y (N) 
      QUOT = 1.0D0 
      S = 0.0D0 
      N1 = 0 
      DO 100 I = 1, N 
         IF (X (I) .NE.0.0D0) THEN 
            S = S + Y (I) / X (I) 
            N1 = N1 + 1 
         ENDIF 
  100 END DO 
      IF (N1.NE.0) THEN 
         QUOT = S / DBLE (N1) 
      ENDIF 
      RETURN 
      END FUNCTION QUOT                             
!                                                                       
!                                                                       
      SUBROUTINE QSCAL (X, Y, R, N) 
!                                                                       
!*****************************************************************      
!*                                                               *      
!*  computes the dot product of two DOUBLE PRECISION vectors     *      
!*                                                               *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (N), Y (N) 
      R = 0.0D0 
      DO 100 I = 1, N 
         R = R + X (I) * Y (I) 
  100 END DO 
      RETURN 
      END SUBROUTINE QSCAL                          
