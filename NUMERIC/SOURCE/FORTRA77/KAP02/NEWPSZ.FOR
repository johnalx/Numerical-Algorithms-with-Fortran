C[BA*)
C[KA{F 2}{Nonlinear Equations in One Variable}
C[       {Nonlinear Equations in One Variable}*)
C[FE{F 2.5.1}
C[  {The Newton Method for Simple Roots}
C[  {The Newton Method for Simple Roots}*)
C[LE*)
      SUBROUTINE NEWPSZ(FCT,FDER,X0,MAXIT,ABSERR,RELERR,X,F,
C[IX{NEWPSZ}*)
     +                  NUMIT,IERR)
C
C*****************************************************************
C                                                                *
C  This SUBROUTINE determines an approximate value X for a zero  *
C  of the function FCT with known derivative FDER by applying    *
C  Newton's method for simple zeros.                             *
C[BE*)
C  The iteration function is:                                    *
C                 PHI(X):=X-FCT(X)/FDER(X).                      *
C  For a simple zero the order of convergence is P=2 with        *
C  efficiency E=1.414...                                         *
C[BA*)
C  If FCT is an algebraic polynomial the                         *
C  SUBROUTINE NEWPOL should be applied instead.                  *
C  For multiple zeros the modified Newton method                 *
C  (SUBROUTINE NEWMOD) is recommended.                           *
C[BE*)
C  Although the order of                                         *
C  convergence is slightly lower for methods that do not use     *
C  derivatives such as SUBROUTINE ZERORF, PEGASU and others,     *
C  they are preferable to Newton's method. Their efficiency is   *
C  generally better since they only require one functional       *
C  evaluation per iteration step.                                *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  FCT      : function for which a zero is to be found.          *
C             It has the form                                    *
C                 DOUBLE PRECISION FUNCTION FCT(X)               *
C             and has to be declared as EXTERNAL within the      *
C             calling program (or as INTRINSIC if a FORTRAN      *
C             standard function is used).                        *
C  FDER     : 1st derivative of the function FCT.                *
C             It is declared as                                  *
C                 DOUBLE PRECISION FUNCTION FDER(X)              *
C             and has to be defined as EXTERNAL within the       *
C             calling program (or as INTRINSIC, if a FORTRAN     *
C             standard function is used).                        *
C  X0       : starting value for the iteration.                  *
C  MAXIT    : maximum number of iteration steps (MAXIT >= 1).    *
C  ABSERR   : ) error parameters. Both have be set to >= 0.0     *
C  RELERR   : ) their sum has to be > 0.0. The following mixed   *
C               test is used:                                    *
C                   ABS(X1-X2) <= ABS(X2)*RELERR+ABSERR.         *
C               Thus if RELERR=0.0, we test for the absolute     *
C               error; if ABSERR=0.0 is chosed, then the         *
C               relative error is tested for.                    *
C               The values entered for ABSERR and RELERR are     *
C               accepted unchanged by the program if they both   *
C               exceed four times the machine constant, or in    *
C               case one of them equals zero, the other exceeds  *
C               four times the machine constant.                 *
C               If this is not the case for both or one of the   *
C               values, they are set to that value internally.   *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  ABSERR   : ) error parameters actually used.                  *
C  RELERR   : )                                                  *
C  X        : last approximate value for the zero of FCT         *
C  F        : function value of FCT at the last approximate      *
C             value X.                                           *
C  NUMIT    : number of iteration steps performed. The number    *
C             of functional evaluations is two times NUMIT.      *
C  IERR     : = 0, input parameter are not as specified.         *
C             = 1, zero found; break-off criterion for the       *
C                  difference between the last two approximations*
C                  has been met.                                 *
C             = 2, zero X found.                                 *
C                      ABS(FCT(X)) < 4 * machine constant.       *
C             = 3, starting value X0 is a zero of FCT:           *
C                      FCT(X0)=0.0 (machine zero).               *
C             = 4, zero not found. The maximum number of         *
C                  iterations was reached without meeting        *
C                  the break-off criterion. Possibly the         *
C                  error parameters were chosen too small or     *
C                  the starting value X0 was not close enough to *
C                  a zero position of FCT.                       *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: MACHPD                                  *
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
C
C  initializing the iteration step counter NUMIT.
C
      NUMIT=0
C
C  testing the validity of the input parameters ABSERR, RELERR and MAXIT.
C
      IF(ABSERR .LT. 0.0D0 .OR. RELERR .LT. 0.0D0 .OR. ABSERR+RELERR
     +   .LE. 0.0D0 .OR. MAXIT .LT. 1) THEN
         IERR=0
         RETURN
      END IF
C
C  test whether X0 already is a zero of FCT:
C        FCT(X0)=0.0 (machine zero).
C
      X=X0
      F=FCT(X)
      IF(F .EQ. 0.0D0) THEN
         IERR=3
         RETURN
      END IF
C
C  computation of the machine constant FMACHP.
C
      FMACHP=1.0D0
   10 FMACHP=0.5D0*FMACHP
      IF(MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
      FMACHP=2.0D0*FMACHP
C
C  in case ABSERR and/or RELERR were chosen too small, they are now set
C  to four times the machine constant.
C
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
C  iteration loop for calculating a new approximate value X for
C  the desired zero. First it we check whether F=FCT(X) is less
C  than four times the machine constant.
C  If this is true the program is ended by setting IERR=2.
C
   20 IF(DABS(F) .LT. DUMMY) THEN
         IERR=2
         RETURN
      ELSE
C
C  calculation of the 1st derivative FD of the function FCT at X;
C  if this turns out to be zero, the derivative will be set to 1.E-6
C  in order to prevent division by zero.
C
         FD=FDER(X)
         IF(DABS(FD) .LT. DUMMY) THEN
            FD=DSIGN(1.0D-6,FD)
         END IF
C
C  increase the iteration step counter by 1 and calculate
C  a new approximation for the zero and its corresponding
C  functional value.
C
         NUMIT=NUMIT+1
         DIFF=F/FD
         X=X-DIFF
         F=FCT(X)
C
C  test whether the break-off criterion is met for the last two
C  approximations or whether the maximum number MAXIT of
C  iterations has been reached.
C
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
