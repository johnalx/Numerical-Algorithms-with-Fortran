C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 4.14}
C[  {Solving Linear Systems via Householder Transformations}
C[  {Solving Linear Systems via Householder Transformations}*)
C[LE*)
      SUBROUTINE SHOUSE (A, LDA, M, N, D, X, IERR)
C[IX{SHOUSE}*)
C[LE*)
C*****************************************************************
C                                                                *
C  SUBROUTINE SHOUSE  solves a linear minimization problem by    *
C  applying Householder transformations, i.e., the euclidean     *
C  norm of A(ORG)*X-B is minimized. Here A is an (M+1)x(N+2)     *
C  matrix A(0:M,0:N+1) with M >= N that contains A(ORG) of       *
C  dimensions (M+1)x(N+1) in the columns 0, ..., N; B(0:M) is a  *
C  vector of length M+1 stored in the (N+1)st column of A, and X *
C  is the solution vector of length N+1.                         *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  A     : 2-dimensional array A(0:LDA,0:N+1) that contains the  *
C          (M+1)x(N+1) matrix A in columns 0 to N, and the       *
C          (M+1)-vector B in column N+1                          *
C  LDA   : leading dimension of A as defined in the calling      *
C          program (LDA has to be >= M)                          *
C  M     : M+1 = number of rows in A (M has to be >= N)          *
C  N     : N+1 = number of columns in A                          *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  X     : (N+1)-vector X(0:N) containing the solution vector    *
C  IERR  : error parameter:                                      *
C          = 0 : everything is o.k.                              *
C          = 1 : numerically the matrix A is not of full rank,   *
C                => no unique solution exists                    *
C                                                                *
C                                                                *
C  AUXILIARY PARAMETER:                                          *
C  ====================                                          *
C                                                                *
C  D     : (N+1)-VECTOR D(0:N), which contains the diagonal      *
C          elements of A during the factorization                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: MACHPD                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Ilona Westermann                                   *
C  date     : 01.09.1987                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LDA, M, N, IERR
      DOUBLE PRECISION A(0:LDA,0:N+1), D(0:N), X(0:N)
      IERR = 0
      FMACHP = 1.0D0
   10 FMACHP = 0.5D0 * FMACHP
      IF (MACHPD(1.0D0 + FMACHP) .EQ. 1) GOTO 10
      FMACHP = 8.0D0 * FMACHP
C
C  Householder transformation:
C
      DO 100 I=0,N
C
C  calculating the essential parameters of the
C  transforming Householder matrices
C
          RADI = A(I,I) * A(I,I)
          DO 40 K=I+1,M
              RADI = RADI + A(K,I) * A(K,I)
   40     CONTINUE
          IF (RADI .LT. FMACHP) THEN
              IERR = 1
              RETURN
          ENDIF
          AIBETR = DSQRT(RADI) * DSIGN(1.0D0,A(I,I))
          AK = 1.0D0 / (RADI + AIBETR * A(I,I))
          A(I,I) = A(I,I) + AIBETR
C
C  update the matrix A and the vector B, stored in the last column of A,
C  by using the new Householder matrix and starting from the left
C
          D(I) = -AIBETR
          DO 100 K=I+1,N+1
              FACTOR = 0.0D0
              DO 50 J=I,M
                  FACTOR = FACTOR + A(J,K) * A(J,I)
   50         CONTINUE
              FACTOR = FACTOR * AK
              DO 100 J=I,M
                  A(J,K) = A(J,K) - FACTOR * A(J,I)
  100 CONTINUE
C
C  determine the solution vector by backsubstitution
C
      DO 80 I=N,0,-1
         SUM = 0.0D0
         DO 70 K=I+1,N
            SUM = SUM + A(I,K) * X(K)
   70    CONTINUE
         X(I) = (A(I,N+1) - SUM) / D(I)
   80 CONTINUE
      RETURN
      END
