C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 3.3.5}{The Laguerre Method}{The Laguerre Method}*)
C[LE*)
      SUBROUTINE LAGUER (A,N,ABSERR,RELERR,MAXIT,XI,NITER,INUM,
C[IX{LAGUER}*)
     +                   WORK,IERR)
C
C*****************************************************************
C                                                                *
C  For a real polynomial PN of degree N, this program computes   *
C  all real roots using Laguerre's method, provided PN has only  *
C  real roots.                                                   *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  A      : (N+1)-vector A(0:N) containing the coefficients      *
C           of the real polynomial PN, where                     *
C              PN(X) = A(0) + A(1)*X + ... + A(N)*X**N           *
C  N      : degree of PN, N > 2                                  *
C  ABSERR : ) error bounds, each of which must be nonnegative,   *
C  RELERR : ) with their sum positive. The following break-off   *
C             criterion is used:                                 *
C                ABS(X1-X2) <= ABS(X2)*RELERR+ABSERR.            *
C             If RELERR=0.0, then we test for the absolute error,*
C             if ABSERR=0.0, then we test the relative error.    *
C             The input values for  ABSERR and RELERR  are used  *
C             without modification only if both exceed four times*
C             the machine constant or , in case one is zero, the *
C             other must exceed that constant. Otherwise on or   *
C             both of the bounds are adjusted internally to that *
C             value.                                             *
C  IMAX   :   maximal number of iterations allowed for each zero,*
C             IMAX >= 1                                          *
C                                                                *
C                                                                *
C  AUXILIARY VARIABLES:                                          *
C  ====================                                          *
C  WORK   : (N+1)-vector WORK(0:N)                               *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  XI     : N-vector XI(1:N) containing the N real roots of PN   *
C  NITER  : INTEGER N-vector NITER(1:N), containing the number of*
C           iterations performed for finding the root with the   *
C           same index                                           *
C  IANZ   : number of roots found                                *
C  IERR   : error parameter                                      *
C           =0, invalid input values for ABSERR, RELERR, IMAX    *
C               or N                                             *
C           =1, all roots have been found                        *
C           =2, maximally allowed number of iterations IMAX has  *
C               been reached; the roots found earlier are stired *
C               in XI(1:IANZ)                                    *
C           =3, the intermediate variable S is negative, i.e.,   *
C               SQRT(S) is not real and there may be complex     *
C               roots of PN                                      *
C           =4, error message from subroutine QUADRA when        *
C               computing the last two roots of PN               *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines: HORN1, HORN2, QUADRA, MACHPD            *
C                                                                *
C*****************************************************************
C                                                                *
C  Author    : Gisela Engeln-MÅllges                             *
C  Date      : 06.03.1992                                        *
C  Source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C  Declarations
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(0:N),WORK(0:N),XI(1:N)
      INTEGER NITER(1:N)
C
C  check input values of  ABSERR, RELERR, IMAX and N
C
      IF(ABSERR .LT. 0.0D0 .OR. RELERR .LT. 0.0D0 .OR.
     +   ABSERR+RELERR .LE. 0.0D0 .OR. MAXIT .LT. 1 .OR.
     +   N .LE. 2) THEN
         IERR=0
         RETURN
      ENDIF
C
C  Compute the machine constant and adjust error bounds if needed
C
      FMACHP=1.0D0
   10 FMACHP=0.5D0*FMACHP
      IF(MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
      FMACHP=2.0D0*FMACHP
      EPS=4.0D0*FMACHP
      IF(RELERR .EQ. 0.0D0) THEN
         IF(ABSERR .LT. EPS) ABSERR=EPS
      ELSEIF(ABSERR .EQ. 0.0D0) THEN
         IF(RELERR .LT. EPS) RELERR=EPS
      ELSE
         IF(ABSERR .LT. EPS) ABSERR=EPS
         IF(RELERR .LT. EPS) RELERR=EPS
      ENDIF
C
C   Initialize
C
      INUM=0
      NACT=N
      IERR=1
C
C  Loop to compute N-2 roots
C
      DO 20 I=1,N-2
C
C  Initialize iteration counter for the I-th root
C
         ITER=0
C
C  Initialize starting value for all roots
C
         X=0.0D0
C
C  Iteration loop
C
   30    IF(ITER .GE. MAXIT) THEN
            IERR=2
            RETURN
         ENDIF
         CALL HORN2(A,N,X,NACT,PN,PN1,PN2,WORK)
         ITER=ITER+1
C
C  Calculate S from the iteration rule and find SQRT(S)
C
         S=(NACT-1)*((NACT-1)*PN1**2-NACT*PN*PN2)
         IF(DABS(S) .LE. EPS) S=0.0D0
         IF(S .LT. 0.0D0) THEN
            IERR=3
            RETURN
         ENDIF
         SRS=DSQRT(S)
C
C  Choose the sign of  SQRT(S) so that it coincides
C  with the sign of PN1
C
         SRS=DSIGN(SRS,PN1)
C
C  Calculate the denominator of the iteration rule.
C  If it is less than four times the machine constant,
C  we set it equal to  EPS
C
         XNENN=PN1+SRS
         IF(DABS(XNENN) .LT. EPS) THEN
            XNENN=DSIGN(EPS,XNENN)
         ENDIF
C
C  Find the difference of successive iterates and a new iterate
C
         DIFF=NACT*PN/XNENN
         X=X-DIFF
C
C  check the break-off criterion
C
         IF(DABS(DIFF) .GT. DABS(X)*RELERR+ABSERR) GOTO 30
         XI(I)=X
         NITER(I)=ITER
         INUM=INUM+1
         CALL HORN1(A,N,X,NACT)
         NACT=NACT-1
   20 CONTINUE
C
C  Finally compute the remaining two roots for the
C  quadratic equation:
C  A(N)*X**2 + A(N-1)*X + A(N-2) = 0
C
      CALL QUADRA(A(N),A(N-1),A(N-2),EPS,X1,X2,IERR1)
      IF(IERR1 .NE. 1) THEN
         IERR=4
         RETURN
      ENDIF
      INUM=INUM+2
      XI(N-1)=X1
      XI(N)=X2
      NITER(N-1)=0
      NITER(N)=0
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE HORN1 (A,N,X0,NACT)
C[IX{HORN1}*)
C
C*****************************************************************
C                                                                *
C  If P has degree NACT, this subroutine finds the remainder     *
C  polynomial PAB of degree NACT-1 that results from dividing off*
C  one known root of P. It is neeeded in the SUBROUTINE LAGUER.  *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  A      : (N+1)-vector A(0:N) containing the coefficients of P:*
C           P(X) = A(N)*X**NACT + A(N-1)*X**(NACT-1) + ...       *
C                               + A(N-NACT+1)*X + A(N-NACT)      *
C  N      : dimension of A + 1, as specified in the calling      *
C           program                                              *
C  X0     : known root of P; P is divided by (X-X0) to obtain PAB*
C  NACT   : actual degree of P                                   *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  A      : (N+1)-vector A(0:N), containing the coefficients of  *
C           the rest polynomial PAB                              *
C           PAB(X) = A(N)*X**(NACT-1) + A(N-1)*X**(NACT-2) + ... *
C                     + A(N-NACT+2)*X + A(N-NACT-1)              *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author    : Gisela Engeln-MÅllges                             *
C  Date      : 06.03.1992                                        *
C  Source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(0:N)
      HELP=0.0D0
      DO 10 K=N,N-NACT,-1
         HELP=HELP*X0+A(K)
         A(K)=HELP
   10 CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE HORN2 (A,N,X,NACT,P,P1,P2,WORK)
C[IX{HORN2}*)
C
C*****************************************************************
C                                                                *
C  This subroutine calculates the function value and that of the *
C  first and second derivative of a given polynom P of degree    *
C  NACT at X. It is needed in SUBROUTINE LAGUER.                 *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  A      : (N+1)-vector A(0:N) containing the coefficients of P:*
C           P(X) = A(N)*X**NACT + A(N-1)*X**(NACT-1) + ...       *
C                               + A(N-NACT+1)*X + A(N-NACT)      *
C  N      : dimension of A + 1, as specified in the calling      *
C           program                                              *
C  X      : value where P(X), P'(X) and P''(X) has to be computed*
C  NACT   : actual degree of P                                   *
C                                                                *
C                                                                *
C  AUXILIARY PARAMETERS:                                         *
C  =====================                                         *
C  WORK   : (N+1)-vector WORK(0:N)                               *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  P      : P(X)                                                 *
C  P1     : P'(X)                                                *
C  P2     : P''(X)                                               *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author    : Gisela Engeln-MÅllges                             *
C  Date      : 02.18.1992                                        *
C  Source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(0:N),WORK(0:N)
      DO 10 I=N-NACT,N
         WORK(I)=A(I)
   10 CONTINUE
      DO 20 K=0,2
         HELP=0.0D0
         DO 30 I=N,N-NACT+K,-1
            HELP=HELP*X+WORK(I)
            WORK(I)=HELP
   30    CONTINUE
   20 CONTINUE
      P=WORK(N-NACT)
      P1=WORK(N-NACT+1)
      P2=2.0D0*WORK(N-NACT+2)
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE QUADRA (A,B,C,EPS,X1,X2,IERR)
C[IX{QUADRA}*)
C
C*****************************************************************
C                                                                *
C  This subroutine solves a real quadratic equation              *
C        A*X**2 + B*X + C = 0  with A not equal to 0,            *
C  that has real roots. It serves the subroutine LAGUER.         *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  A, B, C: real coefficients of the quadratic equation          *
C  EPS    : EPS = 4.0 * machine constant                         *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  X1, X2 : real roots                                           *
C  IERR   : error parameter                                      *
C           =1, all is ok                                        *
C           =2, error: either A=0 or complex roots               *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Authors   : Gisela Engeln-MÅllges                             *
C  Date      : 09.03.1992                                        *
C  Source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IERR=1
      IF(A .EQ. 0.0D0) THEN
         IERR=2
         RETURN
      ENDIF
      D=B*B-4.0D0*A*C
C
C  If the magnitude of the discriminant  D = B**2 - 4*A*C  is less than
C  EPS = 4.0 * machine constant, we set D = 0. If D <  -EPS, we stop:
C  the roots are complex conjugate.
C
      IF(DABS(D) .LE. EPS) D=0.0D0
      IF(D .LT. 0.0D0) THEN
         IERR=2
         RETURN
      ENDIF
      WURZ=DSQRT(D)
C
C  Berechnung der beiden reellen Nullstellen der quadratischen
C  Gleichung
C
      V1=B+WURZ
      V2=B-WURZ
      IF(DABS(V1) .GE. DABS(V2)) THEN
         X1=-2.0D0*C/V1
      ELSE
         X1=-2.0D0*C/V2
      ENDIF
      IF(DABS(X1) .LE. EPS) THEN
         X1=0.0D0
         X2=-(B/A)
      ELSE
         X2=C/(A*X1)
      ENDIF
      RETURN
      END
