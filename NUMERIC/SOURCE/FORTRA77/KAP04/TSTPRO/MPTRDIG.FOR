      PROGRAM  TEST
C                                                ( Thomas Meuser )
C*****************************************************************************
C                                                                            *
C     Testprogram for the subroutines  TRDIG, TRDIGP, TRDIGS.                *
C     Solve a linear system   A*X = RS  with strongly nonsingular tridiagonal*
C     system matrix A.                                                       *
C                                                                            *
C     The test example gives the following output:                           *
C                                                                            *
C[  EXAMPLE:                                                                ]*
C[  ========                                                                ]*
C[  COEFFICIENT MATRIX A:                                                   ]*
C[                                                                          ]*
C[   .50000D+01  -.10000D+01   .00000D+00   .00000D+00                      ]*
C[  -.20000D+01   .50000D+01  -.10000D+01   .00000D+00                      ]*
C[   .00000D+00  -.20000D+01   .30000D+01  -.10000D+01                      ]*
C[   .00000D+00   .00000D+00  -.10000D+01   .50000D+01                      ]*
C[                                                                          ]*
C[  RIGHT HAND SIDE                                                         ]*
C[                                                                          ]*
C[   .40000D+01   .20000D+01   .00000D+00   .40000D+01                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION                                                                ]*
C[                                                                          ]*
C[   .10000D+01   .10000D+01   .10000D+01   .10000D+01                      ]*
C[  STOP. NO ERROR!                                                         ]*
C                                                                            *
C     Other tests possible for different data.                               *
C                                                                            *
C*****************************************************************************
C
      PARAMETER (N=4,M=N*N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION UD(N),HD(N),OD(N),RS(N),X(N),A(N,N)
C
C   Initialize data; change if desired      
C
      DATA (HD(I),I=1,N)/5.D0,5.D0,3.D0,5.D0/
      DATA (UD(I),I=2,N)/-2.D0,-2.D0,-1.D0/
      DATA (OD(I),I=1,N-1)/-1.D0,-1.D0,-1.D0/
      DATA (RS(I),I=1,N)/4.D0,2.D0,0.D0,4.D0/
C
C   Ausgabe des Testbeispiels in Matrixform.
C
      DATA ((A(J,I),I=1,N),J=1,N) /M*0.0D0/
      A(1,1)=HD(1)
      A(1,2)=OD(1)
      DO 100 I=2,N-1
         A(I,I)=HD(I)
         A(I,I+1)=OD(I)
         A(I,I-1)=UD(I)
  100 CONTINUE
      A(N,N)=HD(N)
      A(N,N-1)=UD(N)
      WRITE (*,2010)
      DO 110 I=1,N
  110 WRITE (*,2020)  (A(I,J),J=1,N)
      WRITE (*,2030)  (RS(I),I=1,N)
C
      CALL TRDIG (N,UD,HD,OD,RS,X,MARKE)
C
      WRITE (*,2040)   (X(I),I=1,N)
C
C   Output of error message in  MARKE (-1/0/1)
C
      MARKE=MARKE+2
      GOTO (10,20,30),MARKE
   10 WRITE (*,2050) 'ERROR:  N <= 2 !'
      STOP
   20 WRITE (*,2050)'ERROR: MATRIX A NOT STRONGLY NONSINGULAR!'
      STOP
   30 WRITE (*,2050)'STOP. NO ERROR!'
      STOP
C
 2010 FORMAT (1X,'C[',2X,'EXAMPLE:',T78,']*',/,
     +        1X,'C[',2X,8('='),T78,']*',/,
     +        1X,'C[',2X,'COEFFICIENT MATRIX A:',T78,']*',/,
     +        1X,'C[',T78,']*')
 2020 FORMAT (1X,'C[',4(1X,D12.5),T78,']*')
 2030 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'RIGHT HAND SIDE',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',4(1X,D12.5),T78,']*',/,
     +        1X,'C[',T78,']*')
 2040 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'SOLUTION',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',4(1X,D12.5),T78,']*')
 2050 FORMAT (1X,'C[',2X,A,T78,']*')
      END
