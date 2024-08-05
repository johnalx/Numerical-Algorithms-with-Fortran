![  {Solving Linear Systems via Householder Transformations}            
![  {Solving Linear Systems via Householder Transformations}*)          
      SUBROUTINE SHOUSE (A, LDA, M, N, D, X, IERR) 
!*****************************************************************      
!                                                                *      
!  SUBROUTINE SHOUSE  solves a linear minimization problem by    *      
!  applying Householder transformations, i.e., the euclidean     *      
!  norm of A(ORG)*X-B is minimized. Here A is an (M+1)x(N+2)     *      
!  matrix A(0:M,0:N+1) with M >= N that contains A(ORG) of       *      
!  dimensions (M+1)x(N+1) in the columns 0, ..., N; B(0:M) is a  *      
!  vector of length M+1 stored in the (N+1)st column of A, and X *      
!  is the solution vector of length N+1.                         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A     : 2-dimensional array A(0:LDA,0:N+1) that contains the  *      
!          (M+1)x(N+1) matrix A in columns 0 to N, and the       *      
!          (M+1)-vector B in column N+1                          *      
!  LDA   : leading dimension of A as defined in the calling      *      
!          program (LDA has to be >= M)                          *      
!  M     : M+1 = number of rows in A (M has to be >= N)          *      
!  N     : N+1 = number of columns in A                          *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  X     : (N+1)-vector X(0:N) containing the solution vector    *      
!  IERR  : error parameter:                                      *      
!          = 0 : everything is o.k.                              *      
!          = 1 : numerically the matrix A is not of full rank,   *      
!                => no unique solution exists                    *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETER:                                          *      
!  ====================                                          *      
!                                                                *      
!  D     : (N+1)-VECTOR D(0:N), which contains the diagonal      *      
!          elements of A during the factorization                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Ilona Westermann                                   *      
!  date     : 01.09.1987                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER LDA, M, N, IERR 
      DOUBLEPRECISION A (0:LDA, 0:N + 1), D (0:N), X (0:N) 
      IERR = 0 
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = 8.0D0 * FMACHP 
!                                                                       
!  Householder transformation:                                          
!                                                                       
      DO 100 I = 0, N 
!                                                                       
!  calculating the essential parameters of the                          
!  transforming Householder matrices                                    
!                                                                       
         RADI = A (I, I) * A (I, I) 
         DO 40 K = I + 1, M 
            RADI = RADI + A (K, I) * A (K, I) 
   40    END DO 
         IF (RADI.LT.FMACHP) THEN 
            IERR = 1 
            RETURN 
         ENDIF 
         AIBETR = DSQRT (RADI) * DSIGN (1.0D0, A (I, I) ) 
         AK = 1.0D0 / (RADI + AIBETR * A (I, I) ) 
         A (I, I) = A (I, I) + AIBETR 
!                                                                       
!  update the matrix A and the vector B, stored in the last column of A,
!  by using the new Householder matrix and starting from the left       
!                                                                       
         D (I) = - AIBETR 
         DO 100 K = I + 1, N + 1 
            FACTOR = 0.0D0 
            DO 50 J = I, M 
               FACTOR = FACTOR + A (J, K) * A (J, I) 
   50       END DO 
            FACTOR = FACTOR * AK 
            DO 100 J = I, M 
               A (J, K) = A (J, K) - FACTOR * A (J, I) 
  100 CONTINUE 
!                                                                       
!  determine the solution vector by backsubstitution                    
!                                                                       
      DO 80 I = N, 0, - 1 
         SUM = 0.0D0 
         DO 70 K = I + 1, N 
            SUM = SUM + A (I, K) * X (K) 
   70    END DO 
         X (I) = (A (I, N + 1) - SUM) / D (I) 
   80 END DO 
      RETURN 
      END SUBROUTINE SHOUSE                         
