      PROGRAM TEST
C                                                             89/06/21
C*********************************************************************
C                                                                    *
C  Testprogram one step methods                                      *
C                                                                    *
C  The test example poduces these results:                           *
C[                                                                  ]*
C[ FUNCTION NO :  1     METHOD OF EULER-CAUCHY                      ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .27182807E+01       ]*
C[                                         YEX= .27182818E+01       ]*
C[               IFEHL=   1             ABSERR= .11047761E-05       ]*
C[              NEXT H= .78125000E-03   RELERR= .40642456E-06       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  2 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  2     METHOD OF EULER-CAUCHY                      ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .36787972E+00       ]*
C[                                         YEX= .36787944E+00       ]*
C[               IFEHL=   1             ABSERR= .27510617E-06       ]*
C[              NEXT H= .15625000E-02   RELERR= .74781554E-06       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  2 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  3     METHOD OF EULER-CAUCHY                      ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .19999917E+01       ]*
C[                                         YEX= .20000000E+01       ]*
C[               IFEHL=   1             ABSERR= .82602119E-05       ]*
C[              NEXT H= .19531250E-03   RELERR= .41301230E-05       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  3 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  4     METHOD OF EULER-CAUCHY                      ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .60652997E+00       ]*
C[                                         YEX= .60653066E+00       ]*
C[               IFEHL=   1             ABSERR= .69440437E-06       ]*
C[              NEXT H= .12500000E-01   RELERR= .11448806E-05       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  1 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  5     METHOD OF EULER-CAUCHY                      ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .50000000E+00        Y= .10000000E+01       ]*
C[                                         YEX= .10000000E+01       ]*
C[               IFEHL=   0             ABSERR= .00000000E+00       ]*
C[              NEXT H= .30000000E+00   RELERR= .00000000E+00       ]*
C[ AFTERWARDS :      X= .15000000E+01        Y= .20000000E+01       ]*
C[                                         YEX= .20000000E+01       ]*
C[               IFEHL=   0             ABSERR= .10258461E-12       ]*
C[              NEXT H= .39062500E-03   RELERR= .51292304E-13       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  3 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  1     METHOD OF HEUN                              ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .27182810E+01       ]*
C[                                         YEX= .27182818E+01       ]*
C[               IFEHL=   1             ABSERR= .87167621E-06       ]*
C[              NEXT H= .12500000E-01   RELERR= .32067186E-06       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  2     METHOD OF HEUN                              ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .36787932E+00       ]*
C[                                         YEX= .36787944E+00       ]*
C[               IFEHL=   1             ABSERR= .12156121E-06       ]*
C[              NEXT H= .12500000E-01   RELERR= .33043773E-06       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  3     METHOD OF HEUN                              ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .20000000E+01       ]*
C[                                         YEX= .20000000E+01       ]*
C[               IFEHL=   1             ABSERR= .30390427E-07       ]*
C[              NEXT H= .62500000E-02   RELERR= .15195213E-07       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  4     METHOD OF HEUN                              ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .60653423E+00       ]*
C[                                         YEX= .60653066E+00       ]*
C[               IFEHL=   1             ABSERR= .35665882E-05       ]*
C[              NEXT H= .12500000E-01   RELERR= .58802753E-05       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  5     METHOD OF HEUN                              ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .12500000E+01       ]*
C[                                         YEX= .12500000E+01       ]*
C[               IFEHL=   1             ABSERR= .00000000E+00       ]*
C[              NEXT H= .40000000E+00   RELERR= .00000000E+00       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  1     METHOD OF RUNGE-KUTTA                       ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .27182817E+01       ]*
C[                                         YEX= .27182818E+01       ]*
C[               IFEHL=   1             ABSERR= .17714670E-06       ]*
C[              NEXT H= .10000000E+00   RELERR= .65168633E-07       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  2     METHOD OF RUNGE-KUTTA                       ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .36787941E+00       ]*
C[                                         YEX= .36787944E+00       ]*
C[               IFEHL=   1             ABSERR= .31004156E-07       ]*
C[              NEXT H= .10000000E+00   RELERR= .84278041E-07       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  3     METHOD OF RUNGE-KUTTA                       ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .20000000E+01       ]*
C[                                         YEX= .20000000E+01       ]*
C[               IFEHL=   1             ABSERR= .44408921E-15       ]*
C[              NEXT H= .10000000E+00   RELERR= .22204460E-15       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  4     METHOD OF RUNGE-KUTTA                       ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .60653174E+00       ]*
C[                                         YEX= .60653066E+00       ]*
C[               IFEHL=   0             ABSERR= .10777099E-05       ]*
C[              NEXT H= .20000000E+00   RELERR= .17768400E-05       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C[                                                                  ]*
C[ FUNCTION NO :  5     METHOD OF RUNGE-KUTTA                       ]*
C[ INITIAL VALUES :  X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS :      X= .10000000E+01        Y= .12514286E+01       ]*
C[                                         YEX= .12500000E+01       ]*
C[               IFEHL=   1             ABSERR= .14285714E-02       ]*
C[              NEXT H= .40000000E+00   RELERR= .11415525E-02       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND.                    ]*
C                                                                    *
C*********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FHILF(0:4,3)
      INTEGER IND(3)
      CHARACTER*12 VERF(3)
      LOGICAL JUMP
C
      COMMON /FUNKTI/ IFKT,IMETH,XEND,JUMP
C
      DATA VERF /'EULER-CAUCHY','HEUN        ','RUNGE-KUTTA '/
      JUMP=.TRUE.
C.
C        Schleife Åber die drei Verfahren
C.
      DO 70 IMETH=1,3
C.
C        Intervallende
C.
         XEND = 1.D0
C.
C        LOOP THROUGH THREE TEST METHODS
C.
         DO 50 IFKT=1,5
            CALL BERECH(FHILF,VERF,IND)
50       CONTINUE
70    CONTINUE
      STOP
      END
C
C
      SUBROUTINE BERECH(FHILF,VERF,IND)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FHILF(0:4,3)
      INTEGER IND(3)
      CHARACTER*12 VERF(3)
      LOGICAL JUMP
      COMMON /FUNKTI/ IFKT,IMETH,XEND,JUMP
      EXTERNAL FKT1
      RELEPS = 1.0D-06
      ABSEPS = RELEPS
      H = 0.1D0
      HMAX = 10.D0 * H
      X = 0.D0
      Y = 1.D0
      WRITE (*,3600) IFKT,VERF(IMETH),X,Y,RELEPS
      IFEHL = 0
      IFZAH1 = 0
      IFZAH2 = 0
      IND(1) = 0
      IND(2) = IMETH
      IND(3) = 0
   20 XEND0 = X + 1.D0
C.
C     Test function No 5 has a discontinuous first derivative at 
C     ABS(X)=0.5 ; integrate accordingly
C.
      IF ( IFKT .EQ. 5 .AND. DABS(XEND0) .GT. 0.5D0 .AND. JUMP ) THEN
         XEND0 = 0.5D0
         IND(3) = 1
         JUMP = .NOT.JUMP
      ENDIF
   30 CALL DEQSSP (X,Y,FKT1,XEND0,H,HMAX,ABSEPS,RELEPS,IND,IFEHL,FHILF)
      CALL LFKT1 (X,YEX)
      ABSERR = DABS(Y-YEX)
      RELERR = ABSERR/DABS(Y)
      IF ( IFEHL .LE. 1 ) THEN
         WRITE (*,3700) X, Y, YEX, IFEHL, ABSERR, H, RELERR
C.
C        if function 5 was integrated until the jump in the derivative,
C        integration can proceed normally from there on
C.
         IF ( IFKT .EQ. 5 .AND. DABS(X) .GT. 0.5D0) IND(3) = 0
      ELSEIF ( IFEHL .EQ. 2 ) THEN
         IFZAH1 = IFZAH1 + 1
         GOTO 30
      ELSEIF ( IFEHL .EQ. 3 ) THEN
         H = 0.1D0
         RELEPS = RELEPS*10.D0
         ABSEPS = RELEPS
         IFZAH2 = IFZAH2 + 1
         IFZAH1 = IFZAH1 + 1
         GOTO 30
      ELSE
         WRITE (*,3800) IFEHL
         STOP
      ENDIF
      IF ( X .LT. XEND ) GOTO 20
      WRITE (*,3900) IFZAH1,IFZAH2
      RETURN
 3600 FORMAT (1X,'C[',T70,']*',/,
     +        1X,'C[ FUNCTION NO.:',I3,5X,'METHOD OF ',A,T70,']*',/,
     +        1X,'C[ INITIAL VALUES:   X=',E14.8,8X,'Y=',E14.8,
     +           T70,']*',/,
     +        1X,'C[',38X,'RELEPS=',E14.8,T70,']*')
 3700 FORMAT (1X,'C[ AFTERWARDS:       X=',E14.8,8X,'Y=',E14.8,
     +           T70,']*',/,
     +        1X,'C[',41X,'YEX=',E14.8,T70,']*',/,
     +        1X,'C[',15X,'IFEHL= ',I3,13X,'ABSERR=',E14.8,T70,']*',/,
     +        1X,'C[',14X,'NEXT H=',E14.8,3X,'RELERR=',E14.8,T70,']*')
 3800 FORMAT (1X,'C[ ERRONEOUS INPUT PARAMETER, IFEHL=',I2,T70,']*')
 3900 FORMAT (1X,'C[ AT INTERVAL END:',T70,']*',/,
     +        1X,'C[',15X,I4,' FAILED TRIES',T70,']*',/,
     +        1X,'C[',15X,I4,' REDUCTIONS OF ERROR BOUND.',T70,']*')
      END
C
C    THE TEST FUNCTIONS:
C
      SUBROUTINE FKT1 (X,Y,F)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL JUMP
      COMMON /FUNKTI/ INUM,IMETH,XEND,JUMP
C
      GOTO (10,20,30,40,50) INUM
   10 F = Y
      RETURN
C
   20 F = -Y
      RETURN
C
   30 F = 5.D0*X**4
      RETURN
C
   40 F = -X*Y
      RETURN
C
   50 IF ( DABS(X) .LT. 0.5D0 ) THEN
         F = 0.D0
      ELSEIF ( X .GE. 0.5D0 ) THEN
         F = 2.D0*(X - 0.5D0)
      ELSEIF ( X .LE. -0.5D0 ) THEN
         F = 2.D0*(X + 0.5D0)
      ENDIF
      RETURN
C
      END
C
C    SOLUTIONS OF THE DE's
C
      SUBROUTINE LFKT1 (X,YEX)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL JUMP
      COMMON /FUNKTI/ INUM,IMETH,XEND,JUMP
C
      GOTO (10,20,30,40,50) INUM
C
   10 YEX = DEXP(X)
      RETURN
C
   20 YEX = DEXP(-X)
      RETURN
C
   30 YEX = X**5 + 1.D0
      RETURN
C
   40 YEX = DEXP(-0.5D0*X*X)
      RETURN
C
   50 IF ( DABS(X) .LT. 0.5D0 ) THEN
         YEX = 1.D0
      ELSEIF ( X .GE. 0.5D0 ) THEN
         YEX = (X - 0.5D0)**2 + 1.D0
      ELSEIF ( X .LE. -0.5D0) THEN
         YEX = (X + 0.5D0)**2 + 1.D0
      ENDIF
      RETURN
C
      END
