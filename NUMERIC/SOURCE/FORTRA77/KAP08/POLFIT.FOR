C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 8.1.3.2}
C[  {Discrete Least Squares using Orthogonal Polynomials}
C[  {Discrete Least Squares via Algebraic Polynomials
C[   using Orthogonal Polynomials}*)
C[LE*)
      SUBROUTINE POLFIT (N,X,F,W,NDEG, ALPHA,B,C,SUMERS,IERR)
C[IX{POLFIT}*)
C
C*****************************************************************
C                                                                *
C     This program performs a discrete polynomial approximation  *
C     using orthogonal polynomials.                              *
C[BE*)
C     For N+1 given tuples ( X(J), F(J) ), J = 0, 1, ..., N,     *
C     the polynomial coefficients ALPHA(0), ALPHA(1), .. ,       *
C     ALPHA(NDEG) are computed so that the weighted least squares*
C     error                                                      *
C                                                                *
C      SUMERS =                                                  *
C     ( SUM( J=0 to N ) W(J)*( F(J) - POLY(X(J)) )**2 )** (1/2)  *
C                                                                *
C     becomes minimal for the positive weights W(0), W(1), ..,   *
C     W(N). Here                                                 *
C                                                                *
C       POLY(X) = SUM( K=0 to NDEG ) ALPHA(K) * Q(K,X)           *
C                                                                *
C     denotes the approximating polynomial of degree NDEG.       *
C     The discrete orthogonal polynomials Q(K,X) of degree       *
C     K satisfy the recursion                                    *
C                                                                *
C       Q(K,X) = ( X - B(K) ) * Q(K-1,X) - C(K) * Q(K-2,X),      *
C       for K >= 2 with  Q(0,X) := 1  and  Q(1,X) := X - B(1).   *
C                                                                *
C     By using the values for B(1), B(2), .. , B(NDEG) and       *
C     C(2), .., C(NDEG), the approximating polynomial can be     *
C     evaluated efficiently using the function POEVAL.           *
C     ( the standard FORTRAN conventions apply to the variable   *
C       types used for the transfer parameters).                 *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C        N       : the number of tuples used is N + 1            *
C                  (compare with the dimension of the vectors    *
C                   X, F, W)                                     *
C      X(0:N)    : nodes used                                    *
C      F(0:N)    : functional values at the nodes                *
C      W(0:N)    : weights with W(J) > 0 for all J; for equal    *
C                  weighting set W(J) = 1 for all J.             *
C      NDEG      : degree of the approximating polynomial,       *
C                  0 <= NDEG <= N.                               *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  ALPHA(0:NDEG) : coefficients of the approximating polynomial  *
C                  with respect to the orthogonal polynomials    *
C                  Q(K,X)                                        *
C  B(1:NDEG)   )   coefficients for the recursive calculation    *
C  C(2:NDEG)   )   of the polynomials Q(K,X)                     *
C                  ( if NDEG < 2 or NDEG < 1 the coefficients    *
C                    in C or in B and C are not needed).         *
C  SUMERS        : weighed least squares error of the computed   *
C                  approximating polynomial at the given nodes   *
C  IERR          : error parameter                               *
C                  = 0 : everything o.k.                         *
C                  = 1 : the inequalities 0 <= NDEG <= N are not *
C                        true                                    *
C                  = 2 : not all weights are positive            *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: SCALP, Q, POEVAL                        *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Klaus Niederdrenk                                  *
C  date     : 05.27.87                                           *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(0:N), F(0:N), W(0:N)
      DIMENSION ALPHA(0:NDEG), B(1:NDEG), C(2:NDEG)
C
C ***  checking the input parameters
C
      IERR = 0
      IF ( N .LT. NDEG .OR. NDEG .LT. 0 ) THEN
        IERR = 1
      ELSE
        DO 10 J = 0, N
          IF ( W(J) .LE. 0.0D0 ) IERR = 2
  10    CONTINUE
      ENDIF
      IF ( IERR .NE. 0 ) R E T U R N
C
C *** calculating the coefficients ALPHA(K), B(K) and C(K)
C *** by applying the function SCALP
C
      QKMIN1 = SCALP(0, N, X, F, W, B, C, '(QK,QK)')
      B(1) = SCALP(0, N, X, F, W, B, C, '(XQK,QK)') / QKMIN1
      ALPHA(0) = SCALP(0, N, X, F, W, B, C, '(F,QK)') / QKMIN1
      DO 20 K = 2, NDEG
        QKMIN2 = QKMIN1
        QKMIN1 = SCALP(K-1, N, X, F, W, B, C, '(QK,QK)')
        B(K) = SCALP(K-1, N, X, F, W, B, C, '(XQK,QK)') / QKMIN1
        C(K) = QKMIN1 / QKMIN2
        ALPHA(K-1) = SCALP(K-1, N, X, F, W, B, C, '(F,QK)') / QKMIN1
  20  CONTINUE
      DUMMY = SCALP(NDEG, N, X, F, W, B, C, '(F,QK)')
      ALPHA(NDEG) = DUMMY / SCALP(NDEG, N, X, F, W, B, C, '(QK,QK)')
C
C *** determine the weighted least squares error
C
      SUMERS = 0.0D0
      DO 30 J = 0, N
        DUMMY = F(J) - POEVAL(X(J), NDEG, ALPHA, B, C)
        SUMERS = SUMERS + W(J) * DUMMY * DUMMY
  30  CONTINUE
      SUMERS = DSQRT(SUMERS)
C
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION SCALP (K, N, X, F, W, B, C,CHOICE)
C[IX{SCALP}*)
C
C*****************************************************************
C                                                                *
C     This program uses the formula                              *
C                                                                *
C       ( G1 , G2 ) = SUMME( J=0 to N ) W(J)*G1(X(J))*G2(X(J))   *
C                                                                *
C     to compute the weighted scalar products required for the   *
C     discrete polynomial approximation.                         *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C  subroutines required: Q                                       *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Klaus Niederdrenk                                  *
C  date     : 05.27.87                                           *
C  source   : FORTRAN 77                                         *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(0:N), F(0:N), W(0:N)
      DIMENSION B(1:*), C(2:*)
      CHARACTER *(*) CHOICE
C
      SCALP = 0.0D0
      IF ( CHOICE .EQ. '(QK,QK)' ) THEN
        DO 10 J = 0, N
          DUMMY = Q(K, X(J), B, C)
          SCALP = SCALP + W(J) * DUMMY * DUMMY
  10    CONTINUE
      ELSE IF ( CHOICE .EQ. '(XQK,QK)' ) THEN
        DO 20 J = 0, N
          DUMMY = Q(K, X(J), B, C)
          SCALP = SCALP + W(J) * X(J) * DUMMY * DUMMY
  20    CONTINUE
      ELSE IF ( CHOICE .EQ. '(F,QK)' ) THEN
        DO 30 J = 0, N
          SCALP = SCALP + W(J) * F(J) * Q(K, X(J), B, C)
  30    CONTINUE
      ENDIF
C
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION Q (K, X, B, C)
C[IX{Q}*)
C
C*****************************************************************
C                                                                *
C     This program evaluates the polynomial Q(K,X) of degree K   *
C     at the given node X by applying the two-stage recursion    *
C                                                                *
C       Q(K,X) = ( X - B(K) ) * Q(K-1,X) - C(K) * Q(K-2,X),      *
C       for K >= 2 with Q(0,X) = 1 and Q(1,X) = X - B(1).        *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Klaus Niederdrenk                                  *
C  date     : 05.27.87                                           *
C  source   : FORTRAN 77                                         *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(1:*), C(2:*)
C
      IF ( K .EQ. 0 ) THEN
        Q = 1.0D0
      ELSE IF ( K .EQ. 1 ) THEN
        Q = X - B(1)
      ELSE
        QMIN2 = 1.0D0
        QMIN1 = X - B(1)
        DO 10 I = 2, K
          Q = ( X - B(I) ) * QMIN1 - C(I) * QMIN2
          QMIN2 = QMIN1
          QMIN1 = Q
  10    CONTINUE
      ENDIF
C
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION POEVAL (X, NDEG, ALPHA, B, C)
C[IX{POEVAL}*)
C
C*****************************************************************
C                                                                *
C     In a Horner like way, this program evaluates the           *
C     approximating polynomial                                   *
C                                                                *
C       POLY(X) = SUMME( K=0 to NDEG ) ALPHA(K) * Q(K,X)         *
C                                                                *
C     of degree NDEG, that was determined by SUBROUTINE POLFIT,  *
C     at the given node X.                                       *
C[BE*)
C     ( the standard FORTRAN conventions apply to the variable   *
C       types for the transfer parameters).                      *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C        X       : value at which the approximating polynomial   *
C                  is to be evaluated                            *
C      NDEG      : degree of the approximating polynomial        *
C   ALPHA(0:NDEG): the coefficients of the approximating         *
C                  polynomial determined by SUBROUTINE POLFIT    *
C                  with respect to the orthogonal polynomials    *
C                  Q(K,X)                                        *
C   B(1:NDEG)    : ) coefficients for the recursion              *
C   C(2:NDEG)    : )                                             *
C                                                                *
C  OUTPUT PARAMETER:                                             *
C  =================                                             *
C     POEVAL     : value of the approximating polynomial at X    *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Klaus Niederdrenk                                  *
C  date     : 05.27.87                                           *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ALPHA(0:NDEG), B(1:NDEG), C(2:NDEG)
      IF ( NDEG .EQ. 0 ) THEN
        SK = ALPHA(0)
      ELSE IF ( NDEG .EQ. 1 ) THEN
        SK = ALPHA(0) + ALPHA(1) * ( X - B(1) )
      ELSE
        SKP2 = ALPHA(NDEG)
        SKP1 = ALPHA(NDEG-1) + SKP2 * ( X - B(NDEG) )
        DO 10 K = NDEG-2, 0, -1
          SK = ALPHA(K) + SKP1 * ( X - B(K+1) ) - SKP2 * C(K+2)
          SKP2 = SKP1
          SKP1 = SK
  10    CONTINUE
      ENDIF
      POEVAL = SK
C
      R E T U R N
      END
