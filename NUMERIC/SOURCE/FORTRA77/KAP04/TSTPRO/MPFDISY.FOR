      PROGRAM  TEST
C                                            ( Thomas Meuser )
C*****************************************************************************
C                                                                            *
C     Testprogram for subroutine  FDISY.                                     *
C     Solve a linear system   A*X = RS  with a symmetric positive definite   *
C     five diagonal matrix A.                                                *
C                                                                            *
C     The test example produces the following output:                        *
C                                                                            *
C[  EXAMPLE:                                                                ]*
C[  ========                                                                ]*
C[  COEFFICIENT MATRIX A:                                                   ]*
C[                                                                          ]*
C[   .50000D+01   .10000D+01   .10000D+01   .00000D+00   .00000D+00         ]*
C[   .10000D+01   .50000D+01   .10000D+01   .10000D+01   .00000D+00         ]*
C[   .10000D+01   .10000D+01   .50000D+01   .10000D+01   .10000D+01         ]*
C[   .00000D+00   .10000D+01   .10000D+01   .50000D+01   .10000D+01         ]*
C[   .00000D+00   .00000D+00   .10000D+01   .10000D+01   .50000D+01         ]*
C[                                                                          ]*
C[  RIGHT HAND SIDE                                                         ]*
C[                                                                          ]*
C[   .70000D+01   .80000D+01   .90000D+01   .80000D+01   .70000D+01         ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION                                                                ]*
C[                                                                          ]*
C[   .10000D+01   .10000D+01   .10000D+01   .10000D+01   .10000D+01         ]*
C[  STOP. NO ERROR!                                                         ]*
C                                                                            *
C     Further tests for different test data are possible.                    *
C                                                                            *
C*****************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=5,M=N*N)
      DIMENSION HD(N),OD1(N),OD2(N),RS(N),X(N),A(N,N)
C
C   Initialize test data; change data if desired
C
      DATA (HD(I),I=1,N)/5.D0,5.D0,5.D0,5.D0,5.D0/
      DATA (OD1(I),I=1,N-1)/1.D0,1.D0,1.D0,1.D0/
      DATA (OD2(I),I=1,N-2)/1.D0,1.D0,1.D0/
      DATA (RS(I),I=1,N)/7.D0,8.D0,9.D0,8.D0,7.D0/
C
C   Output in matrix form
C
      DATA ((A(J,I),I=1,N),J=1,N) /M*0.0D0/
      A(1,1)=HD(1)
      A(1,2)=OD1(1)
      A(1,3)=OD2(1)
      A(2,1)=OD1(1)
      A(2,2)=HD(2)
      A(2,3)=OD1(2)
      A(2,4)=OD2(2)
      DO 100 I=3,N-2
         A(I,I)=HD(I)
         A(I,I+1)=OD1(I)
         A(I,I+2)=OD2(I)
         A(I,I-1)=OD1(I-1)
         A(I,I-2)=OD2(I-2)
  100 CONTINUE
      A(N-1,N-1)=HD(N-1)
      A(N-1,N-3)=OD2(N-3)
      A(N-1,N-2)=OD1(N-2)
      A(N-1,N)=OD1(N-1)
      A(N,N)=HD(N)
      A(N,N-1)=OD1(N-1)
      A(N,N-2)=OD2(N-2)
      WRITE (*,2000)
      DO 110 I=1,N
  110 WRITE (*,2010)  (A(I,J),J=1,N)
      WRITE (*,2020)  (RS(I),I=1,N)
C
      CALL FDISY (N,HD,OD1,OD2,RS,X,MARKE)
C
      WRITE (*,2030)   (X(I),I=1,N)
C
C   Output of error code  MARKE (-2/-1/0/1)
C
      MARKE=MARKE+3
      GOTO (10,20,30,40),MARKE
   10 WRITE (*,2040) 'ERROR:  N <= 3 !'
      STOP
   20 WRITE (*,2040)'ERROR: MATRIX A IS NOT POSITIVE DEFINITE!'
      STOP
   30 WRITE (*,2040)'ERROR: MATRIX A NOT STRONGLY NONSINGULAR!'
      STOP
   40 WRITE (*,2040)'STOP. NO ERROR!'
      STOP
C
 2000 FORMAT (1X,'C[',2X,'EXAMPLE:',T78,']*',/,
     +        1X,'C[',2X,8('='),T78,']*',/,
     +        1X,'C[',2X,'COEFFICIENT MATRIX A:',T78,']*',/,
     +        1X,'C[',T78,']*')
 2010 FORMAT (1X,'C[',5(1X,D12.5),T78,']*')
 2020 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'RIGHT HAND SIDE',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',5(1X,D12.5),T78,']*',/,
     +        1X,'C[',T78,']*')
 2030 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'SOLUTION',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',5(1X,D12.5),T78,']*')
 2040 FORMAT (1X,'C[',2X,A,T78,']*')
      END
