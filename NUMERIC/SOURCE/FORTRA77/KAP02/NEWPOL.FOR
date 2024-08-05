C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{}
C[  {The Newton Method for Polynomials}
C[  {The Newton Method for Polynomials}*)
C[LE*)
      SUBROUTINE NEWPOL(A,N,X0,MAXIT,ABSERR,RELERR,X,PN,NUMIT,
C[IX{NEWPOL}*)
     +                  IERR,WORK)
C
C*****************************************************************
C                                                                *
C  This SUBROUTINE determines an approximation X for a zero of   *
C  the algebraic polynomial                                      *
C       PN(X)=A(0)+A(1)*X+A(2)*X**2+...+A(N)*X**N                *
C  with real coefficients A(I) using Newton's method for simple  *
C  roots.                                                        *
C[BE*)
C  The iteration function is:                                    *
C              PHI(X):=X-PN(X)/PNDER(X).                         *
C  Starting with X0, the description of the parameters is        *
C  identical to that of the SUBROUTINE NEWPSZ, as is the organi- *
C  zation of this program with the exception that the functional *
C  evaluations and those for the derivative are determined       *
C  by a Horner scheme (SUBROUTINE NEWPOH).                       *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  A        : (n+1) - vector A(0:N) containing the               *
C             coefficients A(I) of PN.                           *
C  N        : degree of the polynomial PN.                       *
C  X0       : starting value for the iteration.                  *
C  MAXIT    : maximum number of iteration steps (MAXIT >= 1).    *
C  ABSERR   : ) error parameters. Both have to be >= 0.0, while  *
C  RELERR   : ) their sum must be > 0.0. The following mixed     *
C               break-off criterion is used:                     *
C                    ABS(X1-X2) <= ABS(X2)*RELERR+ABSERR.        *
C               Thus if RELERR=0.0 is chosen, this is a test for *
C               the absolute error; if ABSERR=0.0 is chosen, this*
C               is a test for the relative error.                *
C               The values entered for ABSERR and RELERR are     *
C               accepted unchanged by the program if they both   *
C               exceed four times the machine constant, or if one*
C               is zero, if the other exceeds four times the     *
C               machine constant.                                *
C               If this is not the case, both error bounds, or   *
C               one of them are set to that value internally.    *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  ABSERR   : ) error parameters actually used.                  *
C  RELERR   : )                                                  *
C  X        : last approximate value for the zero X              *
C             of PN(X).                                          *
C  PN       : functional value PN(X) at the last approximate X.  *
C  NUMIT    : number of iteration steps executed.                *
C  IERR     : = 0, faulty input parameters.                      *
C             = 1, zero X found; break-off criterion for the     *
C                  difference between the last two approximations*
C                  has been met.                                 *
C             = 2, zero X found.                                 *
C                  ABS(PN(X)) < 4 * machine constant.            *
C             = 3, starting value X0 already is a zero of PN:    *
C                         PN(X0)=0.0 (machine zero).             *
C             = 4, zero has not been found. The maximum number   *
C                  of iterations has been reached without        *
C                  meeting one of the break-off criteria.        *
C                  Possibly the error parameters have been chosen*
C                  too small or the starting value X0 is not     *
C                  close enough to a zero of PN.                 *
C                                                                *
C  AUXILIARY VECTOR:                                             *
C  =================                                             *
C  WORK     : (n+1) - vector WORK(0:N).                          *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: NEWPOH, MACHPD                          *
C                                                                *
C*****************************************************************
C                                                                *
C  author     : Gisela Engeln-Muellges                           *
C  date       : 09.23.1985                                       *
C  source     : FORTRAN 77                                       *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C  declarations.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION A(0:N),WORK(0:N)
C
      NUMIT=0
C
C  test of input parameters.
C
      IF(ABSERR .LT. 0.0D0 .OR. RELERR .LT. 0.0D0 .OR. ABSERR+RELERR
     +   .LE. 0.0D0 .OR. MAXIT .LT. 1) THEN
         IERR=0
         RETURN
      END IF
      X=X0
      CALL NEWPOH(A,N,X,PN,PNDER,WORK)
      IF(PN .EQ. 0.0D0) THEN
         IERR=3
         RETURN
      END IF
      FMACHP=1.0D0
   10 FMACHP=0.5D0*FMACHP
      IF(MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
      FMACHP=2.0D0*FMACHP
      DUMMY=4.0D0*FMACHP
      IF(RELERR .EQ. 0.0D0) THEN
         IF(ABSERR .LT. DUMMY) ABSERR=DUMMY
      ELSE IF(ABSERR .EQ. 0.0D0) THEN
         IF(RELERR .LT. DUMMY) RELERR=DUMMY
      ELSE
         IF(ABSERR .LT. DUMMY) ABSERR=DUMMY
         IF(RELERR .LT. DUMMY) RELERR=DUMMY
      END IF
C
C  iteration loop.
C
   20 IF(DABS(PN) .LT. DUMMY) THEN
         IERR=2
         RETURN
      ELSE
         IF(DABS(PNDER) .LT. DUMMY) THEN
            PNDER=DSIGN(1.0D-6,PNDER)
         END IF
         NUMIT=NUMIT+1
         DIFF=PN/PNDER
         X=X-DIFF
         CALL NEWPOH(A,N,X,PN,PNDER,WORK)
         IF(DABS(DIFF) .LE. DABS(X)*RELERR+ABSERR) THEN
            IERR=1
            RETURN
         ELSE IF(NUMIT .GE. MAXIT) THEN
            IERR=4
            RETURN
         END IF
      END IF
      GOTO 20
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE NEWPOH(A,N,X,PN,PNDER,WORK)
C[IX{NEWPOH}*)
C
c*****************************************************************
C                                                                *
C  This SUBROUTINE uses the Horner scheme in order to calculate  *
C  the functional value PN(X) and the 1st derivative PNDER(X) of *
C  an Nth degree algebraic polynomial with real coefficients     *
C        PN(X)=A(0)+A(1)*X+A(2)*X**2+...+A(N)*X**N               *
C  at X.                                                         *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  A        : (n+1) - vector A(0:N) containing the               *
C             coefficients A(I) of PN.                           *
C  N        : degree of the polynomial PN.                       *
C  X        : value X for which the functional value or the 1st  *
C             derivative value of the polynomial will be         *
C             calculated.                                        *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  PN       : functional value of the polynomial at X.           *
C  PNDER    : 1st derivative of the polynomial at X.             *
C                                                                *
C                                                                *
C  AUXILIARY VECTOR:                                             *
C  =================                                             *
C  WORK     : (n+1) - vector WORK(0:N).                          *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author     : Gisela Engeln-Muellges                           *
C  date       : 09.23.1985                                       *
C  source     : FORTRAN 77                                       *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(0:N),WORK(0:N)
      DO 10 I=0,N,1
         WORK(I)=A(I)
   10 CONTINUE
      DO 20 K=0,1
         S=0.0D0
         DO 20 I=N,K,-1
            S=S*X+WORK(I)
            WORK(I)=S
   20 CONTINUE
      PN=WORK(0)
      PNDER=WORK(1)
      RETURN
      END
