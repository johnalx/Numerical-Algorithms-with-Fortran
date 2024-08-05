      PROGRAM TEST
C
C*****************************************************************************
C                                                                            *
C     Testprogram for the subroutine NEWTIP.                                 *
C                                                                            *
C     Test results from a PC:                                                *
C                                                                            *
C[  TEST EXAMPLE:                                                           ]*
C[  =============                                                           ]*
C[  GIVEN SET OF VALUES:                                                    ]*
C[                                                                          ]*
C[    X(I) I .6000D+01 I .3000D+01 I .0000D+00 I-.3000D+01 I                ]*
C[    -----I-----------I-----------I-----------I-----------I                ]*
C[    Y(I) I  .100D+02 I -.200D+01 I  .200D+01 I -.100D+01 I                ]*
C[                                                                          ]*
C[  OUTPUT OF THE COMPUTED POLYNOMIAL IN THE FOLLOWING FORMAT               ]*
C[                                                                          ]*
C[  N(X) = B(0) + B(1)*(X-X(0))                                             ]*
C[              + B(2)*(X-X(0))*(X-X(1))                                    ]*
C[              ...                                                         ]*
C[              + B(N)*(X-X(0))*...*(X-X(N-1))                              ]*
C[                                                                          ]*
C[    B(  0) =  .100000D+02                                                 ]*
C[    B(  1) =  .400000D+01                                                 ]*
C[    B(  2) =  .888889D+00                                                 ]*
C[    B(  3) =  .141975D+00                                                 ]*
C[                                                                          ]*
C[  FUNCTION VALUE AT  X0= .0000D+00 : F(X)= .2000D+01                      ]*
C[  FUNCTION VALUE AT  X0= .1000D+01 : F(X)= .3086D+00                      ]*
C[  FUNCTION VALUE AT  X0= .2000D+01 : F(X)=-.1309D+01                      ]*
C[  FUNCTION VALUE AT  X0= .3000D+01 : F(X)=-.2000D+01                      ]*
C[  STOP. TERMINATED AS PLANNED!                                            ]*
C                                                                            *
C     Other test runs with differing data are possible.                      *
C                                                                            *
C*****************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (N=3)
      DIMENSION X(0:N),Y(0:N),B(0:N)
      DOUBLE PRECISION NIPFCT
C
C     Initialize test data; change if so desired
C
      DATA (X(I),I=0,N)/6.0D+00,3.0D+00,0.0D+00,-3.0D+00/
      DATA (Y(I),I=0,N)/10.0D+00,-2.0D+00,2.0D+00,-1.0D+00/
C
C     Print out test example
C
      WRITE (*,2000)
      WRITE (*,2010) (X(I),I=0,N)
      WRITE (*,2015)
      WRITE (*,2020) (Y(I),I=0,N)
C
C
      CALL NEWTIP(N,X,Y,B,IFEHL)
C
C
      IF(IFEHL.GT.0)THEN
C
C     Put out error code in  IFEHL (1/2)
C
         GOTO(10,20),IFEHL
   10    WRITE (*,2050)'ERROR: DIMENSION N < 0!'
         STOP
   20    WRITE (*,2050)'THERE ARE TWO IDENTICAL NODES!'
         STOP
      ENDIF
C
C     Output of solution
C
      WRITE (*,2025)
      DO 50 I=0,N
         WRITE (*,2030) I,B(I)
   50 CONTINUE
C
C     Compute functional value at X0
C
      WRITE(*,2035)
      DO 100 I=0,N
         X0=1.0D+00*I
         Y(I)=NIPFCT(X0,X,B,N)
         WRITE(*,2040)X0,Y(I)
  100 CONTINUE
      WRITE (*,2050)'STOP. TERMINATED AS PLANNED!'
      STOP
C
C
 2000 FORMAT (1X,'C[',2X,'TEST EXAMPLE:',T78,']*',/,
     +        1X,'C[',2X,13('='),T78,']*',/,
     +        1X,'C[',2X,'GIVEN SET OF VALUES:',T78,']*',/,
     +        1X,'C[',T78,']*')
 2010 FORMAT (1X,'C[',4X,'X(I)',1X,'I',4(D10.4,1X,'I'),T78,']*')
 2015 FORMAT (1X,'C[',4X,'-----I',4('-----------I'),T78,']*')
 2020 FORMAT (1X,'C[',4X,'Y(I)',1X,'I',4(D10.3,1X,'I'),T78,']*')
 2025 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'OUTPUT OF THE COMPUTED POLYNOMIAL IN ',
     +           'THE FOLLOWING FORMAT',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'N(X) = B(0) + B(1)*(X-X(0))',T78,']*',/,
     +        1X,'C[',14X,'+ B(2)*(X-X(0))*(X-X(1))',T78,']*',/,
     +        1X,'C[',14X,'...',T78,']*',/,
     +        1X,'C[',14X,'+ B(N)*(X-X(0))*...*(X-X(N-1))',T78,']*',/,
     +        1X,'C[',T78,']*')
 2030 FORMAT (1X,'C[',4X,'B(',I3,') = ',D12.6,T78,']*')
 2035 FORMAT (1X,'C[',T78,']*')
 2040 FORMAT (1X,'C[',2X,'FUNCTION VALUE AT  X0=',D10.4,
     +        1X,': F(X)=',D10.4,T78,']*')
 2050 FORMAT (1X,'C[',2X,A,T78,']*')
      END
