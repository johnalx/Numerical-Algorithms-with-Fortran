      PROGRAM  TEST
C
C*****************************************************************************
C                                                                            *
C     Testprogramm for the subroutine CG to solve a linear system A*X = RS   *
C     for a symmetric positive definite  N*N matrix A.                       *
C                                                                            *
C     For the data used, the following output is generated:                  *
C                                                                            *
C[  EXAMPLE:                                                                ]*
C[  ========                                                                ]*
C[  COEFFICIENT MATRIX A:                                                   ]*
C[                                                                          ]*
C[   .50000D+01  -.10000D+01  -.10000D+01  -.10000D+01                      ]*
C[  -.10000D+01   .50000D+01  -.10000D+01  -.10000D+01                      ]*
C[  -.10000D+01  -.10000D+01   .50000D+01  -.10000D+01                      ]*
C[  -.10000D+01  -.10000D+01  -.10000D+01   .50000D+01                      ]*
C[                                                                          ]*
C[  RIGHT HAND SIDE                                                         ]*
C[                                                                          ]*
C[   .20000D+01   .20000D+01   .20000D+01   .20000D+01                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION                                                                ]*
C[                                                                          ]*
C[   .10000D+01   .10000D+01   .10000D+01   .10000D+01                      ]*
C[  STOP. NO ERROR !                                                        ]*
C                                                                            *
C     Further tests with different data are possible.                        *
C                                                                            *
C*****************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (N=4, IA=4, M=N*N)
      PARAMETER (M2=N*(N-1)/2)
      DOUBLE PRECISION A (IA,N), Y (N), D (N), G (N), AMULD (N), X (N)
C
C   initialize; change data for other test examples
C
      DATA (Y(I),I=1,N)/2.0D+00,2.0D+00,2.0D+00,2.0D+00/
C     unteres und oberes Dreieck mit 1, Hauptdiagonale mit 5 vorbesetzen
      DATA ((A(J,I),I=1,J-1),J=2,N)   /M2*-1.0D+00/
      DATA ((A(J,I),I=J+1,N),J=1,N-1) /M2*-1.0D+00/
      DATA (A(I,I),I=1,N) /N*5.0D+00/
C
C   Output of the test matrix
C
      WRITE (*,2000)
      DO 110 I=1,N
        WRITE (*,2010)  (A(I,J),J=1,N)
  110 CONTINUE
      WRITE (*,2020)  (Y(I),I=1,N)
C
      CALL CG (A, N, IA, Y, X, IERR, D, G, AMULD)
C
      WRITE (*,2030)   (X(I),I=1,N)
C
C   Output of error message in IERR (0/1)
C
      IF (IERR .EQ. 0) THEN
        WRITE (*,2040) 'STOP. ERROR --- DIVISION by ZERO in ALPHA!'
        STOP
      ELSE
        WRITE (*,2040) 'STOP. NO ERROR !'
        STOP
      ENDIF
C
 2000 FORMAT (1X,'C[',2X,'EXAMPLE:',T78,']*',/,
     +        1X,'C[',2X,8('='),T78,']*',/,
     +        1X,'C[',2X,'COEFFICIENT MATRIX A:',T78,']*',/,
     +        1X,'C[',T78,']*')
 2010 FORMAT (1X,'C[',4(1X,D12.5),T78,']*')
 2020 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'RIGHT HAND SIDE',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',4(1X,D12.5),T78,']*',/,
     +        1X,'C[',T78,']*')
 2030 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'SOLUTION',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',4(1X,D12.5),T78,']*')
2040  FORMAT (1X,'C[',2X,A,T78,']*')
      END
