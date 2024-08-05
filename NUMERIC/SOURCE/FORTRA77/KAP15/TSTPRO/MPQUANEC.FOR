      PROGRAM TEST
C                                                               93/02/08
C***********************************************************************
C  Test program for SUBROUTINE QUANEC                                  *
C                                                                      *
C                       3.1                                            *
C  Compute the integral  I  3.0*SIN(X) DX   approximately              *
C                       0.0                                            *
C                                                                      *
C                                                                      *
C  Test results:                                                       *
C                                                                      *
C[                                                                    ]*
C[                                                                    ]*
C[ TRAPEZOIDAL RULE :                                                 ]*
C[ LIMITS OF INTEGRATION   A =    0.0000000000D+00                    ]*
C[                         B =    3.1000000000D+00                    ]*
C[                                                                    ]*
C[ METHOD                  1                                          ]*
C[ STEP SIZE               H =    1.2400000000D-01                    ]*
C[ GLOBAL ERROR ORDER :    2                                          ]*
C[                                                                    ]*
C[ VALUE OF INTEGRAL  :    5.9897188053D+00                           ]*
C[                                                                    ]*
C[                                                                    ]*
C[ ERROR ESTIMATE     :    7.6887840528D-03                           ]*
C[                                                                    ]*
C[                                                                    ]*
C[ SIMPSON FORMULA:                                                   ]*
C[ LIMITS OF INTEGRATION   A =    0.0000000000D+00                    ]*
C[                         B =    3.1000000000D+00                    ]*
C[                                                                    ]*
C[ METHOD                  2                                          ]*
C[ STEP SIZE               H =    6.2000000000D-02                    ]*
C[ GLOBAL ERROR ORDER :    4                                          ]*
C[                                                                    ]*
C[ VALUE OF INTEGRAL  :    5.9974059434D+00                           ]*
C[                                                                    ]*
C[                                                                    ]*
C[ ERROR ESTIMATE     :    4.9268429465D-07                           ]*
C[                                                                    ]*
C[                                                                    ]*
C[ 3/8 FORMULA :                                                      ]*
C[ LIMITS OF INTEGRATION   A =    0.0000000000D+00                    ]*
C[                         B =    3.1000000000D+00                    ]*
C[                                                                    ]*
C[ METHOD                  3                                          ]*
C[ STEP SIZE               H =    4.1333333333D-02                    ]*
C[ GLOBAL ERROR ORDER :    4                                          ]*
C[                                                                    ]*
C[ VALUE OF INTEGRAL  :    5.9974056697D+00                           ]*
C[                                                                    ]*
C[                                                                    ]*
C[ ERROR ESTIMATE     :    2.1895334653D-07                           ]*
C[                                                                    ]*
C[                                                                    ]*
C[ 4/90 FORMULA :                                                     ]*
C[ LIMITS OF INTEGRATION   A =    0.0000000000D+00                    ]*
C[                         B =    3.1000000000D+00                    ]*
C[                                                                    ]*
C[ METHOD                  4                                          ]*
C[ STEP SIZE               H =    3.1000000000D-02                    ]*
C[ GLOBAL ERROR ORDER :    6                                          ]*
C[                                                                    ]*
C[ VALUE OF INTEGRAL  :    5.9974054508D+00                           ]*
C[                                                                    ]*
C[                                                                    ]*
C[ ERROR ESTIMATE     :    1.1258466434D-11                           ]*
C[                                                                    ]*
C[                                                                    ]*
C[ 5/288 FORMULA :                                                    ]*
C[ LIMITS OF INTEGRATION   A =    0.0000000000D+00                    ]*
C[                         B =    3.1000000000D+00                    ]*
C[                                                                    ]*
C[ METHOD                  5                                          ]*
C[ STEP SIZE               H =    2.4800000000D-02                    ]*
C[ GLOBAL ERROR ORDER :    6                                          ]*
C[                                                                    ]*
C[ VALUE OF INTEGRAL  :    5.9974054508D+00                           ]*
C[                                                                    ]*
C[                                                                    ]*
C[ ERROR ESTIMATE     :    6.3492886312D-12                           ]*
C[                                                                    ]*
C[                                                                    ]*
C[ 6/840 FORMULA :                                                    ]*
C[ LIMITS OF INTEGRATION   A =    0.0000000000D+00                    ]*
C[                         B =    3.1000000000D+00                    ]*
C[                                                                    ]*
C[ METHOD                  6                                          ]*
C[ STEP SIZE               H =    2.0666666667D-02                    ]*
C[ GLOBAL ERROR ORDER :    8                                          ]*
C[                                                                    ]*
C[ VALUE OF INTEGRAL  :    5.9974054508D+00                           ]*
C[                                                                    ]*
C[                                                                    ]*
C[ ERROR ESTIMATE     :    9.1988326052D-15                           ]*
C[                                                                    ]*
C[                                                                    ]*
C[ 7/17280 FORMULA :                                                  ]*
C[ LIMITS OF INTEGRATION   A =    0.0000000000D+00                    ]*
C[                         B =    3.1000000000D+00                    ]*
C[                                                                    ]*
C[ METHOD                  7                                          ]*
C[ STEP SIZE               H =    1.7714285714D-02                    ]*
C[ GLOBAL ERROR ORDER :    8                                          ]*
C[                                                                    ]*
C[ VALUE OF INTEGRAL  :    5.9974054508D+00                           ]*
C[                                                                    ]*
C[                                                                    ]*
C[ ERROR ESTIMATE     :    2.2997081513D-15                           ]*
C                                                                      *
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER ITEXT*50
      EXTERNAL   FUNK
C
       DO 100 NRV=1,7
         A=0.0D+00
         B=3.1D+00
         N=25
         IF (NRV.EQ.1) THEN
           ITEXT='TRAPEZOIDAL RULE : '
         ELSEIF (NRV.EQ.2) THEN
           ITEXT='SIMPSON FORMULA: '
         ELSEIF (NRV.EQ.3) THEN
           ITEXT='3/8 FORMULA : '
         ELSEIF (NRV.EQ.4) THEN
           ITEXT='4/90 FORMULA : '
         ELSEIF (NRV.EQ.5) THEN
           ITEXT='5/288 FORMULA : '
         ELSEIF (NRV.EQ.6) THEN
           ITEXT='6/840 FORMULA : '
         ELSEIF (NRV.EQ.7) THEN
           ITEXT='7/17280 FORMULA : '
         ENDIF
C
         CALL  QUANEC  (A,B,N,NRV,FUNK,WERT,IFO,FS,IFEHL)
C
         IF (IFEHL .EQ. 0) THEN
           H = (B  -  A) / DBLE(N * NRV)
           WRITE (*,2010) ITEXT,A,B,NRV,H,IFO,WERT
           IF (N.GT.1)  WRITE (*,2020)  FS
         ELSE
           WRITE(*,*) 'ERROR: IFEHL=',IFEHL
           STOP
         ENDIF
100   CONTINUE
      STOP
C
C
 2010 FORMAT (1X, 'C[', T72, ']*', /,
     +        1X, 'C[', T72, ']*', /,
     +        1X, 'C[ ', A, T72, ']*', /,
     +        1X, 'C[ LIMITS OF INTEGRATION', 3X, 'A =', 1PD20.10,
     +            T72, ']*', /,
     +        1X, 'C[', 25X, 'B =', 1PD20.10, T72, ']*', /,
     +        1X, 'C[', T72, ']*', /,
     +        1X, 'C[ METHOD', 14X, I5, T72, ']*', /,
     +        1X, 'C[ STEP SIZE', 15X, 'H =', 1PD20.10,T72, ']*', /,
     +        1X, 'C[ GLOBAL ERROR ORDER :', I5, T72, ']*', /,
     +        1X, 'C[', T72, ']*', /,
     +        1X, 'C[ VALUE OF INTEGRAL  :', 1PD20.10,T72, ']*')
 2020 FORMAT (1X, 'C[', T72, ']*', /,
     +        1X, 'C[', T72, ']*', /,
     +        1X, 'C[ ERROR ESTIMATE     :', 1PD20.10,T72, ']*')
      END
C
      DOUBLE PRECISION FUNCTION FUNK(X)
      DOUBLE PRECISION X
      FUNK = 3.D0*DSIN(X)
      RETURN
      END
