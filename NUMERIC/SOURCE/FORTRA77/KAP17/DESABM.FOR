C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 17.4.3}
C[  {The Predictor-Corrector Method of Adams-Moulton}
C[  {The Predictor-Corrector Method of Adams-Moulton}*)
C[LE*)
      SUBROUTINE DESABM (X,Y,FCT,N,XEND,H,HMAX,ABSERR,RELERR,IND,
     1                   IERR,AUXF,LDF)
C
C*****************************************************************
C                                                                *
C  Numerical solution of a system of ordinary differential       *
C  equations                                                     *
C                                                                *
C      Y' = F(X,Y) with initial condition Y(X0) = Y0             *
C                                                                *
C  using the predictor-corrector method by ADAMS-BASHFORTH-      *
C  MOULTON. The starting values required are produced by the     *
C  RUNGE-KUTTA method with the same order as the A-B-M method.   *
C  An automatic step size control is used, that operates on the  *
C  principle of doubling or halving of the step size with sub-   *
C  sequent error estimation.                                     *
C[BE*)
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C                                                                *
C  X     : starting x value.                                     *
C  Y     : DOUBLE PRECISION vector Y(1:N), the initial values    *
C          at location X.                                        *
C  FCT   : SUBROUTINE, that evaluates the right hand side of the *
C          differential equation. It has to be of the following  *
C          form:                                                 *
C                                                                *
C             SUBROUTINE FCT(X,Y,N,F)                            *
C             DOUBLE PRECISION Y(N), F(N)                        *
C                  .                                             *
C                  .                                             *
C                  .                                             *
C                                                                *
C          In the calling program it must be declared as         *
C          EXTERNAL. Its variables have the following meaning:   *
C             X   independent variable.                          *
C             Y   dependent variable.                            *
C             F   function value of the right hand side of the   *
C                 system of differential equations Y'=F(X,Y).    *
C  N     : number of differential equations in the system.       *
C  XEND  : location where the solution is desired; may be        *
C          smaller than X.                                       *
C  H     : step size for the next step, normally it is           *
C          determined by the program.                            *
C  HMAX  : maximum step size allowed, it has to be positive.     *
C  ABSERR:)                                                      *
C  RELERR:) error parameters, which have to be >= 0. A mixed test*
C           is used:                                             *
C              ABS (LOCAL ERROR) .LE. ABS(Y)*RELERR + ABSERR.    *
C           Thus, if RELERR = 0 is chosen, this corresponds to a *
C           test for the absolute error. If ABSERR = 0 is chosen *
C           this corresponds to a test for the relative error.   *
C           RELERR and ABSERR should be chosen to exceed ten     *
C           times the machine constant. If this is not the case, *
C           they are automatically set equal to this value.      *
C  IND   : indicator, which must be set equal to one on the      *
C          first call.                                           *
C  IERR  : error parameter, on first call it has to be set equal *
C          to  0.                                                *
C  AUXF  : DOUBLE PRECISION auxiliary array AUXF(1:LDF,12),      *
C          LDF .GE. N, where LDF is as defined in the calling    * 
C          program.                                              *
C  LDF   : leading dimension of AUXF.                            *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C                                                                *
C  X     : x value reached during last integration;              *
C          normally X < XEND.                                    *
C  Y     : value for the solution at location X.                 *
C  H     : step size used last, should remain unchanged for the  *
C          next step. If it is changed, IERR has to be reset to  *
C          zero.                                                 *
C  IND   : after the first call set equal to 2, should not be    *
C          changed.                                              *
C  IERR  : error parameter.                                      *
C          = 0: everything o.k., after resetting of XEND         *
C               DESABM can be called again.                      *
C          = 1: everything o.k., after resetting of XEND         *
C               DESABM can be called again.                      *
C               the starting values for the A-B-M-method are     *
C               known, the R-K set up procedure is omitted.      *
C          = 2: after ICTMAX function evaluations the procedure  *
C               did not reach XEND, a repeated call without any  *
C               change of the parameters may be successful.      *
C               (otherwise try to increase the error parameters).*
C               ICTMAX is defined locally as DATA with a preset  *
C               value of 500. Local changes can be made inside   *
C               the program.                                     *
C          = 3: the step size is less than eight times the       *
C               machine constant at an x value. Before further   *
C               calls, H and the error parameters need to be     *
C               increased.                                       *
C          = 4: ABSERR or RELERR is negative or both             *
C               ABSERR = 0.0 and RELERR = 0.0                    *
C          = 5: XEND equals X.                                   *
C          = 6: HMAX is negative.                                *
C                                                                *
C  NOTE:                                                         *
C  =====                                                         *
C                                                                *
C  The point XEND normally is not reached by the program, but due*
C  to the step size control a point X with XEND - H < X < XEND.  *
C  If the solution is required at XEND, DESABM has to be called  *
C  again with H = XEND - X and IERR = 0.                         *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: MACHPD                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Jobst Hoffmann                                     *
C  date     : 07.25.1988 (revised 11.11.1994)                    *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X,XEND,H,HMAX,ABSERR,RELERR
      DIMENSION Y(N),AUXF(LDF,12)
      DIMENSION ALPHA(4),P(4),BETA(4),BETDER(4),FERR(2)
C
      SAVE ALPHA,BETA,BETDER,P,FERR,FMACHP
C
      DATA ICTMAX/500/
C.
C     Initialize the constants. This is performed
C     during the first call of the program only, 
C     afterwards it is omitted.
C.
      IF ( IND .NE. 1 ) GOTO 10
C.
C     machine constant
C.
      FMACHP = 1.0D0
    5 FMACHP = 0.50D0 * FMACHP
      IF (MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 5
      FMACHP = FMACHP * 2.0D0
C.
C     RUNGE-KUTTA coefficients
C.
      ALPHA(1) = 1.0D0/6.0D0
      ALPHA(2) = 1.0D0/3.0D0
      ALPHA(3) = 1.0D0/3.0D0
      ALPHA(4) = 1.0D0/6.0D0
C
      P(1) = 0.0D0
      P(2) = 0.50D0
      P(3) = 0.50D0
      P(4) = 1.0D0
C.
C     ADAMS-BASHFORTH coefficients
C.
      BETA(1) = -9.0D0
      BETA(2) = 37.0D0
      BETA(3) = -59.0D0
      BETA(4) = 55.0D0
C.
C     ADAMS-MOULTON coefficients
C.
      BETDER(1) = 1.0D0
      BETDER(2) = -5.0D0 
      BETDER(3) = 19.0D0 
      BETDER(4) = 9.0D0
C.
C     factors for the error estimate (extrapolation) 
C.
      FERR(1) = 1.0D0/80.0D0
      FERR(2) = -19.0D0/270.0D0 
C 
      IND = 2 
C.
C     Input error check, Initialize variables 
C.
   10 IF ( ABSERR .LT. 0.0 .OR. RELERR .LT. 0.0 .OR.
     1     ABSERR + RELERR .LE. 0.0 ) THEN
         IERR = 4
         RETURN 
      ELSEIF ( XEND .EQ. X ) THEN
         IERR = 5
         RETURN 
      ELSEIF ( HMAX .LE. 0.0D0 ) THEN
         IERR = 6
         RETURN 
      ENDIF 
      IF ( ABSERR .LT. 10.0D0*FMACHP ) ABSERR = 10.0D0*FMACHP 
      IF ( RELERR .LT. 10.0D0*FMACHP ) RELERR = 10.0D0*FMACHP 
      ICOUNT = 0
      IF ( IERR .EQ. 2 ) GOTO 170
      IF ( IERR .EQ. 1 ) GOTO 210
      IERR = 0 
      CALL FCT (X,Y(1),N,AUXF(1,2))
      ICOUNT = ICOUNT + 1 
C.
C     store the starting step size
C.
      HLOC = H
C.
C     determine the current step size
C.
   20 IMETH = 1 
      IF ( HMAX .LT. H ) H = HMAX
      ABSH = DABS(H) 
      H = DSIGN(DMIN1(ABSH,DABS(XEND-X)/3.0D0),XEND - X) 
      ABSH = DABS(H) 
      IF ( ABSH .LE. 8.0D0*FMACHP*DABS(X) ) THEN
         H = HLOC 
         IERR = 0
         RETURN 
      ENDIF 
C.
C     RUNGE-KUTTA steps with error estimates, in order 
C     to determine the starting values for the multi-step 
C     method.
C     If  ISTART = 1:  step size  H. 
C     If  ISTART = 0:  step size 3H. 
C.
   30 ISTART = 1
   40 CONTINUE
C.
C     store the starting values
C.
      X0 = X
      XDUMMY = X 
      DO 50 I=1,N 
         AUXF(I,6) = Y(I)
         AUXF(I,7) = Y(I)
   50 CONTINUE
C.
C     Perform the three steps 
C.
      DO 120 JBEG=1,3 
         DO 60 I=1,N
            AUXF(I,9) = 0.0D0 
   60    CONTINUE 
         IRK = 1
   70    CALL FCT (X0,AUXF(1,7),N,AUXF(1,8))
         DO 80 I=1,N
            AUXF(I,9) = AUXF(I,9) + ALPHA(IRK)*AUXF(I,8) 
   80    CONTINUE 
         IF ( IRK .EQ. 4) GOTO 100
         IRK = IRK + 1
         X0 = XDUMMY + P(IRK)*H
         DO 90 I=1,N
            AUXF(I,7) = AUXF(I,6) + P(IRK)*H*AUXF(I,8) 
   90    CONTINUE 
         GOTO 70
  100    IF ( ISTART .EQ. 0 ) GOTO 130
         XDUMMY = XDUMMY + H
         X0 = XDUMMY 
         DO 110 I=1,N 
            AUXF(I,6) = AUXF(I,6) + H*AUXF(I,9)
            AUXF(I,7) = AUXF(I,6) 
            IF ( JBEG .EQ. 3 ) AUXF(I,11) = AUXF(I,6) 
  110    CONTINUE 
         CALL FCT (X0,AUXF(1,6),N,AUXF(1,JBEG+2)) 
  120 CONTINUE
C.
C     After three steps with step size H, perform one step  
C     with step size 3H 
C.
      ISTART = 1 - ISTART 
      HDUMMY = H 
      H = H + H + H 
      GOTO 40 
  130 ICOUNT = ICOUNT + 13
      DO 140 I=1,N
         AUXF(I,10) = AUXF(I,6) + H*AUXF(I,9)
  140 CONTINUE
      X0 = X + H
C.
C     continue calculations with the original step size
C.
      H = HDUMMY 
  150 IF ( ICOUNT .GE. ICTMAX ) THEN
C.
C        return with message about excessive function calls
C.
         X = X0 
         DO 160 I=1,N 
            Y(I) = AUXF(I,11)
  160    CONTINUE 
         IERR = 2
         RETURN 
      ENDIF 
C.
C     error estimate 
C.
  170 DMAX = 0.0D0 
      YMAX = 0.0D0 
      DO 180 I=1,N
         AUXF(I,12) = FERR(IMETH)*(AUXF(I,11) - AUXF(I,10))
         DMAX = DMAX1(DMAX,DABS(AUXF(I,12)))
         YMAX = DMAX1(YMAX,DABS(AUXF(I,11)))
  180 CONTINUE
      IF ( DMAX .GE. RELERR*YMAX + ABSERR ) THEN
C.
C        unsuccessful run, the step is repeated with half the step size
C.
         H = 0.5D0*H
         IF ( DABS(H) .GT. 8.0D0*FMACHP*DABS(X0) ) GOTO 30 
         IERR = 3
         RETURN 
      ENDIF 
C.
C     the step was successful, calculation is continued with the
C     old step size
C.
      X = X0
      DO 190 I=1,N
         AUXF(I,11) = AUXF(I,11) + AUXF(I,12)
         Y(I) = AUXF(I,11) 
  190 CONTINUE
      CALL FCT (X0,Y(1),N,AUXF(1,5)) 
      ICOUNT = ICOUNT + 1 
      IF ( DMAX .LE. 0.02D0*( RELERR*YMAX + ABSERR ) ) THEN 
C.
C        The precision obtained is too high, the step size is 
C        doubled, continue with the RUNGE-KUTTA starting procedure
C.
         H = 2.0D0*H 
         HLOC = DMAX1(H,HLOC) 
         DO 200 I=1,N 
            AUXF(I,2) = AUXF(I,5) 
  200    CONTINUE 
         IF ( H .GT. 0. .AND. X0 .LT. XEND .OR. 
     1        H .LT. 0. .AND. X0 .GT. XEND ) GOTO 20
         IERR = 0
         RETURN 
      ENDIF 
      IF ( H .GT. 0. .AND. X0 + H .GE. XEND .OR.
     1     H .LT. 0. .AND. X0 + H .LE. XEND ) THEN
C.
C        We are at the end of the interval. For further calls IERR is set to 
C        1, so that on repeated calls the calculations can be continued
C        with ADAMS-BASHFORTH-MOULTON steps
C.
         IF ( ABSH .LE. 8.0D0*FMACHP*DABS(X0) ) THEN
            H = HLOC
            IERR = 0 
            RETURN
         ENDIF
         IERR = 1
         RETURN 
      ENDIF 
C.
C     start of the ADAMS-BASHFORTH-MOULTON steps
C.
  210 IMETH = 2 
      DO 230 J=1,4
         DO 220 I=1,N 
            AUXF(I,J) = AUXF(I,J+1) 
  220    CONTINUE 
  230 CONTINUE
C.
C     predictor-step 
C.
      DO 240 I=1,N
         AUXF(I,8) = 0.0D0
  240 CONTINUE
      DO 260 J=1,4
         DO 250 I=1,N
            AUXF(I,8) = AUXF(I,8) + BETA(J)*AUXF(I,J)
  250    CONTINUE 
  260 CONTINUE
      DO 270 I=1,N
         AUXF(I,10) = AUXF(I,11) + H*AUXF(I,8)/24.0D0 
  270 CONTINUE
C.
C     new node
C.
      X0 = X0 + H 
      CALL FCT (X0,AUXF(1,10),N,AUXF(1,5))
C. 
C     corrector step
C.
      DO 280 I=1,N
         AUXF(I,8) = 0.0D0
  280 CONTINUE
      DO 300 J=1,4
         DO 290 I=1,N 
            AUXF(I,8) = AUXF(I,8) + BETDER(J)*AUXF(I,J+1)
  290    CONTINUE 
  300 CONTINUE
      DO 310 I=1,N
         AUXF(I,11) = AUXF(I,11) + H*AUXF(I,8)/24.0D0 
  310 CONTINUE
C.
C     function value of the corrector 
C.
      CALL FCT (X0,AUXF(1,11),N,AUXF(1,5))
      ICOUNT = ICOUNT + 2 
      GOTO 150
      END 
