C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 2.8.5}{Zeroin Method}{Zeroin Method}*)
C[LE*)
      SUBROUTINE ZEROIN(FCT,ABSERR,RELERR,MAXFCT,ITAPE,A,B,FB,
C[IX{ZEROIN}*)
     +                  NUMFCT,IERR)
C
C*****************************************************************
C                                                                *
C  The SUBROUTINE ZEROIN finds a zero of odd order of a          *
C  continuous function FCT in the interval [A,B] provided that   *
C  FCT(A) and FCT(B) have different signs, i.e.,                 *
C  FCT(A)*FCT(B) < 0.                                            *
C  The Zeroin method combines the bisection and the secant method*
C  with inverse quadratic interpolation.                         *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  FCT    : Function, whose zero one wants to find. It has the   *
C           form   DOUBLE PRECISION FUNCTION FCT(X)  and must be *
C           declared as EXTERNAL in the calling program, or as   *
C           INTRINSIC, if it is described by standard FORTRAN-77 *
C           functions.                                           *
C  ABSERR : ) error bounds with ABSERR >= 0 and RELERR >= 0.     *
C  RELERR : ) Their sum must be positive. The break-off criterion*
C             used is as follows:                                *
C             ABS(XM) <= 0.5*(ABSERR+ABS(B)*RELERR), where       *
C             XM denotes half the interval length, XM = (B-A)/2. *
C             If RELERR=0, then we only test the absolute error; *
C             if ABSERR=0, we only test the relative arror.      *
C             The input values for ABSERR and RELERR are only    *
C             used by the subroutine if each exceeds four times  *
C             the machine constant, or if one is zero then the   *
C             other must exceed that bound. Otherwise one or both*
C             of them are internally adjusted to four times the  *
C             machine constant.                                  *
C  MAXFCT : Maximal number of function evaluations allowed       *
C  ITAPE  : > 0, Number for a data set that will absorb inter-   *
C                mediate results                                 *
C  A, B   : endpoints of the interval that contains a zero of FCT*
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  ABSERR : ) actually used error bounds                         *
C  RELERR : )                                                    *
C  B      : approximate zero                                     *
C  FB     : functionmal value at the approximate zero B          *
C  NUMFCT : number of actual functional evaluations performed    *
C  IERR   : error parameter:                                     *
C           =-2, ABSERR or RELERR is negative or both are zero   *
C                or  MAXFCT < 1                                  *
C           =-1, FCT(A)*FCT(B) < 0.0 is not true                 *
C           = 0, A or B are numerical zeros of FCT               *
C           = 1, B is a zero with FCT(B)=0.0                     *
C           = 2, the desired accuracy has been achieved:         *
C                ABS(XM) <= error bound                          *
C           = 3, maximal number of function evaluations has been *
C                reached without meeting the break-off criterion *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines: MACHPD, BI                              *
C                                                                *
C*****************************************************************
C                                                                *
C  Authors     : Siegmar Neuner, Gisela Engeln-M�llges           *
C  Date        : 06.01.1992                                      *
C  Source      : FORTRAN 77                                      *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      DOUBLE PRECISION FMACHP,FA,FB,FC,A,B,C,D,E,ABSERR,RELERR,
     +                 TOL1,EPS,XM,R,Q,S,P,FCT,HELP
      INTEGER MAXFCT,ITAPE,NUMFCT,IERR,MACHPD

C
C  Compute four times the machine constant FMACHP
C
      FMACHP=1.0D0
   10 FMACHP=0.5D0*FMACHP
      IF(MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
      FMACHP=2.0D0*FMACHP
      EPS=4.0D0*FMACHP
C
C  Compute FCT at A, B
C
      FA=FCT(A)
      FB=FCT(B)
      NUMFCT=2
C
C  Check whether  FCT(A)*FCT(B) < 0.0D0
C
      HELP=FA*FB
      IF(HELP .GT. 0.0D0) THEN
         IERR=-1
         RETURN
      ELSE IF(HELP .EQ. 0.0D0) THEN
         IERR=0
         RETURN
      END IF
C
C  Check input eror parameters
C
      IF(ABSERR .LT. 0.0D0 .OR. RELERR .LT. 0.0D0 .OR.
     +   ABSERR+RELERR .LE. 0.0D0 .OR. MAXFCT .LT. 1) THEN
         IERR=-2
         RETURN
      END IF
      IF(RELERR .EQ. 0.0D0) THEN
         IF(ABSERR .LT. EPS) ABSERR=EPS
      ELSE IF(ABSERR .EQ. 0.0D0) THEN
         IF(RELERR .LT. EPS) RELERR=EPS
      ELSE
         IF(ABSERR .LT. EPS) ABSERR=EPS
         IF(RELERR .LT. EPS) RELERR=EPS
      END IF
C
C  No zero between  B and C, set C so that there is a
C  zero between  B and C
C
   20 C=A
      FC=FA
      D=B-A
      E=D
C
C  If  FC is the smaller sized function value,
C  swap the interval ends
C
   30 IF(DABS(FC) .LT. DABS(FB)) THEN
         A=B
         B=C
         C=A
         FA=FB
         FB=FC
         FC=FA
      END IF
C
C  TOL1 is an auxiliary variable used for the mixed error test
C
      TOL1=0.5D0*(ABSERR+RELERR*DABS(B))
C
C  take half the interval length XM
C
      XM=0.5D0*(C-B)
      IF(ITAPE .GT. 0) THEN
         WRITE(ITAPE,900) A,B,C
         WRITE(ITAPE,910) FA,FB,FC
      END IF
      R=0.0D0
C
C  If XM is less than TOL1, we have achieved a sufficiently
C  good approximate zero
C
      IF(DABS(XM) .LE. TOL1) THEN
         IERR=2
         RETURN
C
C  Check whether the value  FB  of the best approximation to
C  a zero already vanishes
C
      ELSE IF(FB .EQ. 0.0D0) THEN
         IERR=1
         RETURN
      END IF
      IF(DABS(E) .LT. TOL1) THEN
         CALL BI(XM,ITAPE,D,E)
      ELSE
         IF(DABS(FA) .LE. DABS(FB)) THEN
            CALL BI(XM,ITAPE,D,E)
         ELSE
C
C  If  A and C are different, then together with B one can use these three
C  points for an inverse quadratse interpolation
C
            IF(A .NE. C) THEN
               Q=FA/FC
               R=FB/FC
               S=FB/FA
               P=S*(2.0D0*XM*Q*(Q-R)-(B-A)*(R-1.0D0))
               Q=(Q-1.0D0)*(R-1.0D0)*(S-1.0D0)
            ELSE
C
C  Here we use the secant method or interpolate linearly
C
               S=FB/FA
               P=2.0D0*XM*S
               Q=1.0D0-S
            END IF
C
C  The sign of P/Q is reversed for the following division
C
            IF(P .GT. 0.0D0) THEN
               Q=-Q
            ELSE
               P=DABS(P)
            END IF
            IF((2.0D0*P) .GE. (3.0D0*XM*Q-DABS(TOL1*Q))) THEN
               CALL BI(XM,ITAPE,D,E)
            ELSE
               IF(P .GE. DABS(0.5D0*E*Q)) THEN
                  CALL BI(XM,ITAPE,D,E)
               ELSE
C
C  For either interpolation we compute the quotient P/Q
C  which shall be added to B
C
                  E=D
                  D=P/Q
                  IF(ITAPE .GT. 0) THEN
                     IF(R .EQ. 0.0D0) THEN
                        WRITE(ITAPE,*) 'Secant method'
                     ELSE
                        WRITE(ITAPE,*) 'Inverse quadratic ',
     +                                 'interpolation'
                     END IF
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C  The previous best approximate zero B is stored in A
C  and the function value FB in FA
C
      A=B
      FA=FB
C
C  If D exceeds TOL1, it is added to B
C
      IF(DABS(D) .GT. TOL1) THEN
         IF(ITAPE .GT. 0) WRITE(ITAPE,920) D
         B=B+D
      ELSE
C
C  The desired accuracy has been achieved.
C  The best approximate zero B is improved by adding the error bound.
C
         IF(ITAPE .GT. 0) THEN
            IF(XM .LT. 0.0D0) THEN
               WRITE(ITAPE,930) TOL1
            ELSE
               WRITE(ITAPE,940) TOL1
            END IF
         END IF
         B=B+DSIGN(TOL1,XM)
      END IF
C
C  Compute the new value FB at B
C
      FB=FCT(B)
C
C  The iterationcounter is upped by 1
C
      NUMFCT=NUMFCT+1
      IF(ITAPE .GT. 0) WRITE(ITAPE,950) B,FB,NUMFCT
      IF(NUMFCT .GT. MAXFCT) THEN
         IERR=3
         RETURN
      END IF
C
C  If the signs of the function at  B and C are opposite, then
C  there is a zero between B and C
C
      IF((FB*(FC/DABS(FC))) .GT. 0.0D0) GOTO 20
      GOTO 30
C
  900 FORMAT(1X,'A = ',D20.14,'  B = ',D20.14,'  C = ',D20.14)
  910 FORMAT(1X,'FA= ',D20.14,'  FB= ',D20.14,'  FC= ',D20.14)
  920 FORMAT(1X,'distance to the new B:  D= ',D20.14)
  930 FORMAT(1X,'error bound is subtracted:  D= -',D20.14)
  940 FORMAT(1X,'error bound is added:  D= ',D20.14)
  950 FORMAT(1X,'B = ',D20.14,'  FB= ',D20.14,/,
     +       1X,'number of function evaluations: ',I4,//)
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BI(XM,ITAPE,D,E)
C[IX{BI}*)
C
C*****************************************************************
C                                                                *
C  The SUBROUTINE BI is a bisection algorithm for the SUBROUTINE *
C  ZEROIN.                                                       *
C                                                                *
C*****************************************************************
C[BE*)
C
      DOUBLE PRECISION XM,D,E
      D=XM
      E=D
      IF(ITAPE .GT. 0) WRITE(ITAPE,*) 'Bisection'
      RETURN
      END
