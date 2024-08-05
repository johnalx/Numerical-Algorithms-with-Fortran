      PROGRAM TEST
C                                                            4/21/90
C*********************************************************************
C                                                                    *
C  Testprogram for the  Adams-Bashforth-Moulton method               *
C                                                                    *
C  The following data is computed:                                   *
C                                                                    *
C[                                                                  ]*
C[ FUNCTION NO.:  1                                                 ]*
C[ INITIAL VALUES:   X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS:       X= .10000000E+01        Y= .27182818E+01       ]*
C[                                         YEX= .27182818E+01       ]*
C[               IFEHL=   0             ABSERR= .74916359E-07       ]*
C[              NEXT H= .20093879E+02   RELERR= .27560189E-07       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND                     ]*
C[                                                                  ]*
C[ FUNCTION NO.:  2                                                 ]*
C[ INITIAL VALUES:   X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS:       X= .10000000E+01        Y= .36787947E+00       ]*
C[                                         YEX= .36787944E+00       ]*
C[               IFEHL=   0             ABSERR= .28929918E-07       ]*
C[              NEXT H= .20093879E+02   RELERR= .78639663E-07       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND                     ]*
C[                                                                  ]*
C[ FUNCTION NO.:  3                                                 ]*
C[ INITIAL VALUES:   X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS:       X= .10000000E+01        Y= .20000000E+01       ]*
C[                                         YEX= .20000000E+01       ]*
C[               IFEHL=   0             ABSERR= .26041669E-07       ]*
C[              NEXT H= .48225309E+01   RELERR= .13020834E-07       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND                     ]*
C[                                                                  ]*
C[ FUNCTION NO.:  4                                                 ]*
C[ INITIAL VALUES:   X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS:       X= .10000000E+01        Y= .60652976E+00       ]*
C[                                         YEX= .60653066E+00       ]*
C[               IFEHL=   0             ABSERR= .89591077E-06       ]*
C[              NEXT H= .48225309E+01   RELERR= .14771093E-05       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND                     ]*
C[                                                                  ]*
C[ FUNCTION NO.:  5                                                 ]*
C[ INITIAL VALUES:   X= .00000000E+00        Y= .10000000E+01       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS:       X= .10000000E+01        Y= .12500000E+01       ]*
C[                                         YEX= .12500000E+01       ]*
C[               IFEHL=   0             ABSERR= .00000000E+00       ]*
C[              NEXT H= .66979595E+02   RELERR= .00000000E+00       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND                     ]*
C[                                                                  ]*
C[ FUNCTION NO.:  6                                                 ]*
C[ INITIAL VALUES:   X= .00000000E+00     Y(1)= .10000000E+01       ]*
C[                                        Y(2)= .00000000E+00       ]*
C[                                      RELEPS= .10000000E-05       ]*
C[ AFTERWARDS:       X= .10000000E+01     Y(1)= .54030233E+00       ]*
C[                                        Y(2)= .84147095E+00       ]*
C[                                      YEX(1)= .54030231E+00       ]*
C[                                      YEX(2)= .84147098E+00       ]*
C[               IFEHL=   0             ABSERR= .39755802E-07       ]*
C[              NEXT H= .20093879E+02   RELERR= .47245603E-07       ]*
C[ AT INTERVAL END:                                                 ]*
C[                  0 FAILED TRIES                                  ]*
C[                  0 REDUCTIONS OF ERROR BOUND                     ]*
C                                                                    *
C*********************************************************************
C
      PARAMETER (NDIM=2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(NDIM),FHILF(NDIM,12),YEX(NDIM)
      LOGICAL JUMP
C
      COMMON /FUNKTI/ IFKT,XEND,JUMP
C
      JUMP=.TRUE.
C
C.
C     Interval end
C.
      XEND = 1.D0
C.
C     loop over test functions
C.
      DO 60 IFKT=1,6
         CALL BERECH(NDIM,Y,FHILF,YEX)
60    CONTINUE
      STOP
      END
C
C
      SUBROUTINE BERECH(NDIM,Y,FHILF,YEX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(NDIM),FHILF(NDIM,12),YEX(NDIM)
      LOGICAL JUMP
      COMMON /FUNKTI/ IFKT,XEND,JUMP
      EXTERNAL FKTN
      RELEPS = 1.D-06
      ABSEPS = RELEPS
      H = 0.1D0
      HMAX = 10.D0 * H
      X = 0.D0
      Y(1) = 1.D0
      Y(2) = 0.D0
      N = 1
      IF ( IFKT .EQ. 6 ) N = 2
      IF ( N .EQ. 1 ) THEN
         WRITE (*,4100) IFKT,X,Y(1),RELEPS
      ELSEIF ( N .EQ. 2 ) THEN
         WRITE (*,4200) IFKT,X,Y(1),Y(2),RELEPS
      ENDIF
      IFEHL = 0
      IFZAH1 = 0
      IFZAH2 = 0
   20 XEND0 = X + 1.D0
C.
C     Test function No 5 has a discontinuous derivative at ABS(X)=0.5,
C     integrate accordingly
C.
      IF ( IFKT .EQ. 5 .AND. DABS(XEND0) .GT. 0.5D0 .AND. JUMP ) THEN
         H = 0.5D0 - X
         IFEHL = 0
         JUMP = .NOT.JUMP
      ENDIF
   30 CALL DESEXT (X,Y,FKTN,N,XEND0,H,HMAX,ABSEPS,RELEPS,
     F             IFEHL,FHILF,NDIM)
      IF ( IFEHL .LE. 1 ) THEN
         CALL LFKTN (X,YEX,N)
         ABSERR = 0.D0
         DO 40 I=1,N
            ABSERR = DMAX1(ABSERR,DABS(Y(I)-YEX(I)))
            RELERR = ABSERR/DABS(Y(I))
   40    CONTINUE
         IF ( N .EQ. 1 ) THEN
            WRITE (*,4300) X,Y(1),YEX(1),IFEHL,ABSERR,H,RELERR
         ELSEIF ( N .EQ. 2 ) THEN
            WRITE (*,4400) X,Y(1),Y(2),YEX(1),YEX(2),IFEHL,ABSERR,
     F                     H,RELERR
         ENDIF
      ELSEIF ( IFEHL .EQ. 2 ) THEN
         WRITE (*,4500) IFEHL
         RETURN
      ELSE
         WRITE (*,3800) IFEHL
         RETURN
      ENDIF
      IF ( X .LT. XEND ) GOTO 20
      WRITE (*,3900) IFZAH1,IFZAH2
      RETURN
 4100 FORMAT (1X,'C[',T70,']*',/,
     +        1X,'C[ FUNCTION NO.:',I3,T70,']*',/,
     +        1X,'C[ INITIAL VALUES:   X=',E14.8,8X,'Y=',E14.8,
     +           T70,']*',/,
     +        1X,'C[',38X,'RELEPS=',E14.8,T70,']*')
 4200 FORMAT (1X,'C[',T70,']*',/,
     +        1X,'C[ FUNCTION NO.:',I3,T70,']*',/,
     +        1X,'C[ INITIAL VALUES:   X=',E14.8,5X,'Y(1)=',E14.8,
     +           T70,']*',/,
     +        1X,'C[',40X,'Y(2)=',E14.8,T70,']*',/,
     +        1X,'C[',38X,'RELEPS=',E14.8,T70,']*')
 4300 FORMAT (1X,'C[ AFTERWARDS:       X=',E14.8,8X,'Y=',E14.8,
     +           T70,']*',/,
     +        1X,'C[',41X,'YEX=',E14.8,T70,']*',/,
     +        1X,'C[',15X,'IFEHL= ',I3,13X,'ABSERR=',E14.8,T70,']*',/,
     +        1X,'C[',14X,'NEXT H=',E14.8,3X,'RELERR=',E14.8,T70,']*')
 4400 FORMAT (1X,'C[ AFTERWARDS:       X=',E14.8,5X,'Y(1)=',E14.8,
     +           T70,']*',/,
     +        1X,'C[',40X,'Y(2)=',E14.8,T70,']*',/,
     +        1X,'C[',38X,'YEX(1)=',E14.8,T70,']*',/,
     +        1X,'C[',38X,'YEX(2)=',E14.8,T70,']*',/,
     +        1X,'C[',15X,'IFEHL= ',I3,13X,'ABSERR=',E14.8,T70,']*',/,
     +        1X,'C[',14X,'NEXT H=',E14.8,3X,'RELERR=',E14.8,T70,']*')
 3800 FORMAT (1X,'C[ ERRONEOUS INPUT PARAMETER, IFEHL=',I2,T70,']*')
 3900 FORMAT (1X,'C[ AT INTERVAL END:',T70,']*',/,
     +        1X,'C[',15X,I4,' FAILED TRIES',T70,']*',/,
     +        1X,'C[',15X,I4,' REDUCTIONS OF ERROR BOUND',T70,']*')
 4500 FORMAT (1X,'C[ STEPSIZE LESS THAN FOUR TIMES MACHINE',T70,']*',/,
     +        1X,'C[ CONSTANT, IFEHL=',I2,T70,']*')
      END
C
C
C    VARIOUS RIGHT HAND SIDES FOR THE DE's, N-DIMENSIONAL SYSTEMS
C
      SUBROUTINE FKTN (X,Y,N,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),F(N)
      LOGICAL JUMP
C
      COMMON /FUNKTI/ INUM,XEND,JUMP
C
      GOTO (10,20,30,40,50,60) INUM
   10 F(1) = Y(1)
      RETURN
C
   20 F(1) = -Y(1)
      RETURN
C
   30 F(1) = 5.D0*X**4
      RETURN
C
   40 F(1) = -X*Y(1)
      RETURN
C
   50 IF ( DABS(X) .LT. 0.5D0 ) THEN
         F(1) = 0.D0
      ELSEIF ( X .GE. 0.5D0 ) THEN
         F(1) = 2.D0*(X - 0.5D0)
      ELSEIF ( X .LE. -0.5D0 ) THEN
         F(1) = 2.D0*(X + 0.5D0)
      ENDIF
      RETURN
C
   60 F(1) = -Y(2)
      F(2) =  Y(1)
      RETURN
C
      END
C
C
C   SOLUTIONS FOR THE DE's
C
      SUBROUTINE LFKTN (X,YEX,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YEX(N)
      LOGICAL JUMP
C
      COMMON /FUNKTI/ INUM,XEND,JUMP
C
      GOTO (10,20,30,40,50,60) INUM
C
   10 YEX(1) = DEXP(X)
      RETURN
C
   20 YEX(1) = DEXP(-X)
      RETURN
C
   30 YEX(1) = X**5 + 1.D0
      RETURN
C
   40 YEX(1) = DEXP(-0.5D0*X*X)
      RETURN
C
   50 IF ( DABS(X) .LT. 0.5D0 ) THEN
         YEX(1) = 1.D0
      ELSEIF ( X .GE. 0.5D0 ) THEN
         YEX(1) = (X - 0.5D0)**2 + 1.D0
      ELSEIF ( X .LE. -0.5D0) THEN
         YEX(1) = (X + 0.5D0)**2 + 1.D0
      ENDIF
      RETURN
C
   60 YEX(1) = DCOS(X)
      YEX(2) = DSIN(X)
      RETURN
C
      END
