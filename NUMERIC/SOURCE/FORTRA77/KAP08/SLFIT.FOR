C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 8.1.3.4}
C[  {Householder Transformation for Linear Least Squares}
C[  {Solving Linear Least Squares Problems using Householder
C[   Transformations}*)
C[LE*)
      SUBROUTINE SLFIT (X, Y, W, IWFL, FCT, PSI, LDA, M, N, INUM,
C[IX{SLFIT}*)
     +                  A, D, C, FV, SQERR, IERR)
C
C*****************************************************************
C                                                                *
C  SUBROUTINE SLFIT solves the linear discrete least squares     *
C  problem for M+1 nodes (X(I), Y(I)) by computing a linear      *
C  approximation function                                        *
C         F(X) = C(0) * F1(X) + ... + C(N) * FN(X).              *
C  The model functions F1(X), ..., FN(X) must be provided by the *
C  user in form of a SUBROUTINE.                                 *
C[BE*)
C  SLFIT determines the optimal coefficients C(I), I=0, ..., N   *
C  in the least squares sense.                                   *
C  The problem is reduced to a linear minimization problem which *
C  is solved by applying Householder transformations to the      *
C  system matrix of the normal equations.                        *
C  For INUM user-given values the approximating function is then *
C  evaluated.                                                    *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C                                                                *
C  X      (N+1)-vector X(0:M) containing the X-values of the     *
C         nodes                                                  *
C  Y      (N+1)-vector Y(0:M) containing the Y-values at the     *
C         nodes                                                  *
C  W      (N+1)-vector W(0:M) containing the positive weights    *
C  IWFL   if IWFL = 0, the nodes are weighed according to W      *
C         Otherwise the nodes are all equally weighed by 1, and  *
C         W does not need to be defined explicitly.              *
C  FCT    A SUBROUTINE that has to be provided by the user. It   *
C         defines the model functions and has the form:          *
C                                                                *
C         SUBROUTINE FCT (X, N, F)                               *
C         IMPLICIT DOUBLE PRECISION (A-H,O-Z)                    *
C         INTEGER N                                              *
C         DIMENSION F(0:N), X                                    *
C         -----------------                                      *
C         F(0) = value for the 0-th model function at X          *
C         .                                                      *
C         .                                                      *
C         .                                                      *
C         F(N) = value of the N-th model function at X           *
C         -----------------                                      *
C         RETURN                                                 *
C         END                                                    *
C                                                                *
C         In the calling program FCT has to be declared as       *
C         EXTERNAL.                                              *
C  PSI    INUM-vector PSI(1:INUM) containing the values where    *
C         the approximating function is to be evaluated.         *
C  LDA    leading dimension of array A as defined in the calling *
C         program; LDA has to be >= M.                           *
C  M      M+1 = number of nodes                                  *
C  N      N+1 = number of model functions                        *
C  INUM   for INUM values the approximation function is to be    *
C         evaluated.                                             *
C                                                                *
C                                                                *
C  AUXILIARY PARAMETERS:                                         *
C  =====================                                         *
C                                                                *
C  A      2-dim. array A(0:LDA,0:N+1) for SUBROUTINE SHOUSE      *
C  D      (N+1)-vector D(0:N) for SUBROUTINE SHOUSE              *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C                                                                *
C  C      (N+1)-vector C(0:N) containing the optimal coefficients*
C         with respect to the least square error                 *
C  FV     INUM-vector FV(1:INUM) containing the values  of the   *
C         approximating function at the INUM given locations PSI.*
C  IERR  error parameter                                         *
C         = 0  everything  o.k.                                  *
C         = 1  error in the input parameters                     *
C         = 2  approximating function cannot be determined,      *
C              since numerically A does not have maximal rank.   *
C              I.e., the model functions are numerically linearly*
C              dependent.                                        *
C  SQERR least square error                                      *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: SLPRE, SHOUSE, SENORM                   *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Ilona Westermann                                   *
C  date     : 09.01.1987                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(0:M), Y(0:M), A(0:LDA,0:N+1), D(0:N), C(0:N),
     +     PSI(INUM), FV(INUM), W(0:M)
      EXTERNAL FCT
      IERR = 0
C
C     checking the input parameters
C
      IF (M .GE. N .AND. LDA .GE. M .AND. N .GE. 0 .AND. M .GE. 0)
     +   THEN
C
C        minimize (A * X - B) via the following steps:
C
C        1. determine the matrix A in SUBROUTINE SLPRE
C
         CALL SLPRE (X, W, IWFL, FCT, LDA, M, N, D, A)
C
C        2. store the vector B (if necessary modified by the
C           weights) in the last column of A
C
         IF (IWFL .EQ. 0) THEN
            DO 20 I=0,M
               A(I,N+1) = Y(I) * DSQRT(W(I))
   20       CONTINUE
         ELSE
            DO 30 I=0,M
               A(I,N+1) = Y(I)
   30       CONTINUE
         ENDIF
C
C        3. solve the minimization problem using the SUBROUTINE SHOUSE
C
         CALL SHOUSE (A, LDA, M, N, D, C, MARK)
C
C        test for singularity
C
         IF (MARK .EQ. 0) THEN

C
C           determine the value of the approximating function
C           at the locations desired
C
            DO 40 I=1,INUM
               CALL FCT (PSI(I), N, D)
               FV(I) = 0.0D0
               DO 40 J=0,N
                  FV(I) = FV(I) + C(J) * D(J)
   40       CONTINUE
C
C           determine least square error
C
            SQERR = DSQRT (SENORM (A(N+1,N+1), M-(N+1)))
         ELSE
            IERR = 2
         ENDIF
      ELSE
         IERR = 1
      ENDIF
      RETURN
      END
