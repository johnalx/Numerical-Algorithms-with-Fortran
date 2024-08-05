      PROGRAM TEST
C                                                ( Thomas Meuser )
C************************************************************************
C                                                                       *
C     Testprogram for the subroutines  CYCTR, CYCTRP, CYCTRS.           *
C     Solve a linear system  A*X = RS  with cyclic tridiagonal strongly *
C     nonsingular matrix A.                                             *
C                                                                       *
C     The test data produces the following output:                      *
C                                                                       *
C[  EXAMPLE:                                                           ]*
C[  ========                                                           ]*
C[  COEFFICIENT MATRIX A:                                              ]*
C[                                                                     ]*
C[   .50000D+01  -.10000D+01   .00000D+00  -.10000D+01                 ]*
C[   .30000D+01   .20000D+01   .50000D+01   .00000D+00                 ]*
C[   .00000D+00   .70000D+01   .30000D+01  -.10000D+01                 ]*
C[  -.30000D+01   .00000D+00   .80000D+01   .20000D+01                 ]*
C[                                                                     ]*
C[  RIGHT HAND SIDE                                                    ]*
C[                                                                     ]*
C[   .30000D+01   .10000D+02   .90000D+01   .70000D+01                 ]*
C[                                                                     ]*
C[  SOLUTION                                                           ]*
C[                                                                     ]*
C[   .10000D+01   .10000D+01   .10000D+01   .10000D+01                 ]*
C[  STOP. NO ERROR!                                                    ]*
C                                                                       *
C     Other tests are possible with other data.                         *
C                                                                       *
C************************************************************************
C
      PARAMETER(N=4,M=N*N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION UD(N),HD(N),OD(N),UZ(N),RSP(N),RS(N),X(N),A(N,N)
C
C   Initialize data; start further tests with different data
C
      DATA(OD(I),I=1,N-1)/-1.D0,5.D0,-1.D0/
      DATA(HD(I),I=1,N)/5.D0,2.D0,3.D0,2.D0/
      DATA(UD(I),I=2,N)/3.D0,7.D0,8.D0/
      DATA UZ(1)/-3.D0/,RSP(1)/-1.D0/
      DATA(RS(I),I=1,N)/3.D0,10.D0,9.D0,7.D0/
C
C   Output of test example in matrix form
C
      DATA((A(J,I),I=1,N),J=1,N)/M*0.0D0/
      A(1,1)=HD(1)
      A(1,2)=OD(1)
      DO 110 I=2,N-1
          A(I,I)=HD(I)
          A(I,I+1)=OD(I)
          A(I,I-1)=UD(I)
  110 CONTINUE
      A(N,N)=HD(N)
      A(N,N-1)=UD(N)
      A(1,N)=RSP(1)
      A(N,1)=UZ(1)
      WRITE(*,2000)
      DO 120 I=1,N
  120 WRITE(*,2010) (A(I,J),J=1,N)
      WRITE(*,2020) (RS(I),I=1,N)
C
      CALL CYCTR(N,UD,HD,OD,UZ,RSP,RS,X,MARKE)
C
      WRITE(*,2030) (X(I),I=1,N)
C
C   Output of  MARKE (-1,0,1)
C
      MARKE=MARKE+2
      GOTO (10,20,30),MARKE
   10 WRITE(*,2040) 'ERROR:  N <= 2 !'
      STOP
   20 WRITE(*,2040) 'ERROR: MATRIX A NOT STRONGLY NONSINGULAR!'
      STOP
   30 WRITE(*,2040) 'STOP. NO ERROR!'
      STOP
C
 2000 FORMAT (1X,'C[',2X,'EXAMPLE:',T73,']*',/,
     +        1X,'C[',2X,8('='),T73,']*',/,
     +        1X,'C[',2X,'COEFFICIENT MATRIX A:',T73,']*',/,
     +        1X,'C[',T73,']*')
 2010 FORMAT (1X,'C[',4(1X,D12.5),T73,']*')
 2020 FORMAT (1X,'C[',T73,']*',/,
     +        1X,'C[',2X,'RIGHT HAND SIDE',T73,']*',/,
     +        1X,'C[',T73,']*',/,
     +        1X,'C[',4(1X,D12.5),T73,']*')
 2030 FORMAT (1X,'C[',T73,']*',/,
     +        1X,'C[',2X,'SOLUTION',T73,']*',/,
     +        1X,'C[',T73,']*',/,
     +        1X,'C[',4(1X,D12.5),T73,']*')
 2040 FORMAT (1X,'C[',2X,A,T73,']*')
      END
