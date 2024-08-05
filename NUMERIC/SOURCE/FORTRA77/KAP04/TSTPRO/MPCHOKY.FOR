      PROGRAM  TEST
C                                                ( Thomas Meuser )
C************************************************************************
C                                                                       *
C     Testprogramm for the subroutine CHOKY to solve a linear syste     *
C     A*X = RS  for a symmetric positive definite  N*N matrix A via     *
C     the root free Cholesky method.                                    *
C                                                                       *
C                                                                       *
C[  EXAMPLE:                                                           ]*
C[  ========                                                           ]*
C[  COEFFICIENT MATRIX A:                                              ]*
C[                                                                     ]*
C[   .50000D+01  -.10000D+01  -.10000D+01  -.10000D+01                 ]*
C[  -.10000D+01   .50000D+01  -.10000D+01  -.10000D+01                 ]*
C[  -.10000D+01  -.10000D+01   .50000D+01  -.10000D+01                 ]*
C[  -.10000D+01  -.10000D+01  -.10000D+01   .50000D+01                 ]*
C[                                                                     ]*
C[  RIGHT HAND SIDE                                                    ]*
C[                                                                     ]*
C[   .20000D+01   .20000D+01   .20000D+01   .20000D+01                 ]*
C[                                                                     ]*
C[  SOLUTION                                                           ]*
C[                                                                     ]*
C[   .10000D+01   .10000D+01   .10000D+01   .10000D+01                 ]*
C[  STOP. NO ERROR!                                                    ]*
C                                                                       *
C     Further tests with different data are possible.                   *
C                                                                       *
C************************************************************************
C
      PARAMETER (N=4,IA=N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  A(IA,N),X(N),RS(N),Z(N)
C
C   Initialize data; different test examples for different input data
C
      DATA ((A(J,I),I=1,N),J=1,N) /5.D0,-1.D0,-1.D0,-1.D0,-1.D0,5.D0
     +,-1.D0,-1.D0,-1.D0,-1.D0,5.D0,-1.D0,-1.D0,-1.D0,-1.D0,5.D0/
      DATA (RS(I),I=1,N) /2.D0,2.D0,2.D0,2.D0/
C
C   Output of the test example
C
      WRITE (*,2010)
      DO 100  I=1,N
  100 WRITE (*,2020)  (A(I,J),J=1,N)
      WRITE (*,2030)  (RS(I),I=1,N)
C
      CALL CHOKY (N,A,IA,RS,X,Z,MARKE)
C
      WRITE (*,2040)  (X(I),I=1,N)
C
C   Output of error message according to  MARKE (-1/0/1)
C
      MARKE=MARKE+2
      GOTO (10,20,30),MARKE
   10 WRITE(*,2050) 'ERROR: MATRIX A IS NOT POSITIVE DEFINITE!'
      STOP
   20 WRITE(*,2050) 'ERROR: MATRIX A IS NOT STRONGLY NONSINGULAR!'
      STOP
   30 WRITE(*,2050) 'STOP. NO ERROR!'
      STOP
C
 2010 FORMAT (1X,'C[',2X,'EXAMPLE:',T73,']*',/,
     +        1X,'C[',2X,8('='),T73,']*',/,
     +        1X,'C[',2X,'COEFFICIENT MATRIX A:',T73,']*',/,
     +        1X,'C[',T73,']*')
 2020 FORMAT (1X,'C[',4(1X,D12.5),T73,']*')
 2030 FORMAT (1X,'C[',T73,']*',/,
     +        1X,'C[',2X,'RIGHT HAND SIDE',T73,']*',/,
     +        1X,'C[',T73,']*',/,
     +        1X,'C[',4(1X,D12.5),T73,']*')
 2040 FORMAT (1X,'C[',T73,']*',/,
     +        1X,'C[',2X,'SOLUTION',T73,']*',/,
     +        1X,'C[',T73,']*',/,
     +        1X,'C[',4(1X,D12.5),T73,']*')
 2050 FORMAT (1X,'C[',2X,A,T73,']*')
      END
