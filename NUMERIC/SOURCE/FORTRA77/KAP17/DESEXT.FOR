C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 17.5}
C[  {Bulirsch-Stoer-Gragg Extrapolation}
C[  {Bulirsch-Stoer-Gragg Extrapolation}*)
C[LE*)
      SUBROUTINE DESEXT (X,Y,FCT,N,XEND,H,HMAX,ABSERR,RELERR,
C[IX{DESEXT}*)
     F                   IERR,AUXF,LDF)
C
C*****************************************************************
C                                                                *
C  Numerical solution of a system of ordinary differential       *
C  equations                                                     *
C                                                                *
C      Y' = F(X,Y)  with  given initial condition Y(X0) = Y0     *
C                                                                *
C  using the extrapolations method of BULIRSCH-STOER.            *
C  The maximal order of the extrapolation is determined depending*
C  the machine constant; the step size control lollows [HALL76]. *
C[BE*)
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C                                                                *
C  X     : initial x value.                                      *
C  Y     : REAL vector Y(1:N), the initial value at X.           *
C  FCT   : SUBROUTINE, that evaluates the right hand side of the *
C          system of differential equations. It has to be of the *
C          following form:                                       *
C                                                                *
C             SUBROUTINE FCT(X,Y,F,N)                            *
C             DOUBLE PRECISION X, Y(N), F(N)                     *
C                  .                                             *
C                  .                                             *
C                  .                                             *
C                                                                *
C          In the calling program it has to be defined as        *
C          EXTERNAL. In it the variables have the following      *
C          meaning:                                              *
C             X   independent variable.                          *
C             Y   dependent variable.                            *
C             F   function value of the right hand side of the   *
C                 system of differential equations Y'=F(X,Y)     *
C             N   number of differential equations.              *
C  N     : number of differential equations.                     *
C  XEND  : location where the solution is desired; may be        *
C          smaller than X.                                       *
C  H     : step size for the next step. Normally it is deter-    *
C          mined by the program.                                 *
C  HMAX  : maximally allowed step size, it has to be positivs.   *
C  ABSERR:)                                                      *
C  RELERR:) error parameters, which have to be >= 0. A mixed     *
C           test is performed:                                   *
C              ABS (LOCAL ERROR) .LE. ABS(Y)*RELERR + ABSERR.    *
C           Thus if RELERR = 0 is chosen, this corresponds to a  *
C           test for the absolute error. If ABSERR = 0 is chosen,*
C           this corresponds to a test for the relative error.   *
C           RELERR and ABSERR should be chosen larger than ten   *
C           times the machine constant. If this is not the case, *
C           they are automatically set equal to this value.      *
C  IERR  : error parameter, on first call it has to be set equal *
C          to 0.                                                 *
C  AUXF  : REAL auxiliary array AUXF(1:LDF,12), LDF .GE. N, with *
C          LDF as defined in the calling program.                *
C  LDF   : leading dimension of AUXF.                            *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C                                                                *
C  X     : x value that was reached after last during integration*
C          Normally  X = XEND.                                   *
C  Y     : solution value at X.                                  *
C  H     : step size used last, this should remain unchanged for *
C          the next step. If it is changed IERR has to be set to *
C          equal zero.                                           *
C  IERR  : error parameter.                                      *
C          = 0: everything o.k., after resetting of XEND,        *
C               DESEXT can be called again.                      *
C          = 1: after 700 function evaluations the procedure did *
C               not reach XEND. A repeated call without any      *
C               change of the parameters may be successful.      *
C               (otherwise try with increased error parameters)  *
C          = 2: the step size is below four times the machine    *
C               constant for an x value. Before further calls H  *
C               and the error parameters must be increased.      *
C          = 3: ABSERR or RELERR is negative or both             *
C               ABSERR = 0.0 and RELERR = 0.0.                   *
C          = 4: XEND is equal to X.                              *
C          = 5: HMAX is negative.                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: MACHPD                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Jobst Hoffmann                                     *
C  date     : 04.21.1990                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION Y(N),AUXF(LDF,12),BUFOL(12)
      LOGICAL IND
C
      SAVE IEXTMX, FMACHP, IND
C
      DATA IND /.TRUE./, ICTMAX/700/
      DATA BUFOL/ 2.0D0, 4.0D0, 6.0D0, 8.0D0,12.0D0, 16.0D0,
     F           24.0D0,32.0D0,48.0D0,64.0D0,96.0D0,128.0D0/
      IF ( IND ) THEN
C.
C     determine machine constant
C.
        FMACHP = 1.0D0
    5   FMACHP = 0.5D0 * FMACHP
        IF (MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 5
        FMACHP = FMACHP * 2.0D0
C.
C     determine the maximal level of extrapolation. Due to
C     increasing rounding errors effects, extrapolation should not be
C     performed at a higher level.
C.
        IEXTMX = INT(-DLOG(FMACHP)/DLOG(2.0D0)/7.0D0 + 0.5D0)
        IND = .FALSE.
      ENDIF
C.
C     Checking input data, initializing variables.
C.
      IF ( ABSERR .LT. 0.0D0 .OR. RELERR .LT. 0.0D0 .OR.
     F     ABSERR + RELERR .LE. 0.0D0 ) THEN
         IERR = 3
         RETURN
      ELSEIF ( XEND .EQ. X ) THEN
         IERR = 4
         RETURN
      ELSEIF ( HMAX .LE. 0.0D0 ) THEN
         IERR = 5
         RETURN
      ENDIF
      ICOUNT = 0
      IF ( IERR .EQ. 1 ) GOTO 150
      ILINE = 0
      YMAX = 0.0D0
C.
C     determin the first step size
C.
   20 ABSH = DABS(H)
      IF ( HMAX .LT. ABSH ) ABSH = HMAX
      HLOC = DSIGN(DMIN1(ABSH,DABS(XEND-X)),XEND - X)
      IF ( DABS(HLOC) .LE. FMACHP*DABS(X) ) THEN
         IERR = 0
         RETURN
      ENDIF
C.
C     determine the step size for the extrapolation.
C.
   30 ILINE = ILINE + 1
      H0 = HLOC/BUFOL(ILINE)
C.
C     EULER step, store the initial values.
C.
      X0 = X
      DO 40 I=1,N
         AUXF(I,9) = Y(I)
   40 CONTINUE
      CALL FCT (X0,AUXF(1,9),N,AUXF(1,12))
      DO 50 I=1,N
         AUXF(I,10) = AUXF(I,9) + H0*AUXF(I,12)
   50 CONTINUE
      X0 = X0 + H0
      CALL FCT (X0,AUXF(1,10),N,AUXF(1,12))
C.
C     the midpoint rule is applied
C     BUFOL(ILINE)-1 - times.
C.
      IBUFO = INT(BUFOL(ILINE)) - 1
      DO 90 NLOOPS=1,IBUFO
         DO 60 I=1,N
            AUXF(I,11) = AUXF(I,9) + 2.0D0*H0*AUXF(I,12)
   60    CONTINUE
         X0 = X0 + H0
         CALL FCT (X0,AUXF(1,11),N,AUXF(1,12))
C.
C     alter storage for the next step.
C.
         DO 80 J=9,11
            DO 70 I=1,N
               AUXF(I,J) = AUXF(I,J+1)
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
C.
C     stabilize using the trapezoidal rule.
C.
      CALL FCT (X0,AUXF(1,10),N,AUXF(1,12))
      DO 100 I=1,N
         AUXF(I,ILINE) = 0.5D0*(AUXF(I,10) + AUXF(I,9) +
     F                         H0*AUXF(I,12))
  100 CONTINUE
C.
C     at least two values are required for the extrapolation.
C.
      IF ( ILINE .EQ. 1 ) GOTO 30
C.
C     extrapolation.
C.
      MINZ = MIN0(ILINE,IEXTMX)
      DO 120 ICOLUM=2,MINZ
         IDUMMY = MIN0(11-ICOLUM,ILINE-ICOLUM+2)
         DMAX = 0.0D0
         DO 110 I=1,N
            AUXF(I,IDUMMY-1) = AUXF(I,IDUMMY) + ( AUXF(I,IDUMMY)
     F       - AUXF(I,IDUMMY-1))/((BUFOL(ILINE)/BUFOL(IDUMMY-1))**2 - 1)
            YMAX = DMAX1(YMAX,DABS(AUXF(I,IDUMMY-1)))
            DMAX = DMAX1(DMAX,DABS(AUXF(I,IDUMMY-1) - AUXF(I,IDUMMY)))
  110    CONTINUE
         IF ( DMAX .LT. RELERR*YMAX + ABSERR ) GOTO 180
  120 CONTINUE
      ICOUNT = ICOUNT + INT(BUFOL(ILINE)) + 1
      IF ( ICOUNT .GE. ICTMAX ) THEN
C.
C        return with message for excessive functional evaluations,
C        store the current data.
C.
         IERR = 1
         AUXF(1,10) = X
         DO 130 I=1,N
            AUXF(I,9) = Y(I)
  130    CONTINUE
         X = X0
         DO 140 I=1,N
            Y(I) = AUXF(I,IDUMMY-1)
  140    CONTINUE
         RETURN
      ENDIF
      GOTO 170
C.
C     entry point on repeat call following return due to
C     excessive functional evaluations.
C.
  150    X = AUXF(1,10)
         DO 160 I=1,N
            Y(I) = AUXF(I,9)
  160    CONTINUE
  170 IF ( ILINE .LT. 8 ) GOTO 30
C.
C     despite a completed extrapolation scheme the
C     required accuracyprecision was not achieved,
C     calculations are repeated with a smaller step size.
C.
      H = 0.9D0*0.6D0**(ILINE-IEXTMX)*H
      IF ( DABS(H) .LE. 4.0D0*FMACHP*DABS(X0) ) THEN
         IERR = 2
         RETURN
      ENDIF
      ILINE = 0
      GOTO 30
C.
C     the required accuracy was achieved.
C.
  180 X = X + HLOC
      DO 190 I=1,N
         Y(I) = AUXF(I,IDUMMY-1)
  190 CONTINUE
C.
C     the next step size is determined.
C.
      H = 0.9D0*0.6D0**(ICOLUM - IEXTMX)*H
      ILINE = 0
      GOTO 20
      END
