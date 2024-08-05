      PROGRAM TEST
C                                             (Thomas Meuser)     2/11/87
C*************************************************************************
C                                                                        *
C     Test program for the subroutines  ORTOGP, ORTPOL                   *
C     We compute the nodes and weights of generalized Gaussian quadrature*
C                                                                        *
C     The test example produces these results:                           *
C                                                                        *
C[                                                                      ]*
C[ TEST EXAMPLE:                                                        ]*
C[ =============                                                        ]*
C[                                                                      ]*
C[ GIVEN INTEGRAL VALUES:            I        AINT(I)                   ]*
C[                               -------------------------              ]*
C[                                   0  I         1.00000               ]*
C[                                   1  I         1.00000               ]*
C[                                   2  I         2.00000               ]*
C[                                   3  I         6.00000               ]*
C[                                   4  I        24.00000               ]*
C[                                   5  I       120.00000               ]*
C[                                                                      ]*
C[ RESULTS:                                                             ]*
C[     I        NODES X(I)        WEIGHTS (I)                           ]*
C[ --------------------------------------------                         ]*
C[     1  I     .4157745568  I     .7110930099                          ]*
C[     2  I    2.2942803603  I     .2785177336                          ]*
C[     3  I    6.2899450829  I     .0103892565                          ]*
C[                                                                      ]*
C[ NO ERROR                                                             ]*
C[                                                                      ]*
C                                                                        *
C     Other tests with other data are possible.                          *
C                                                                        *
C     Remark: The given values of the integral are derived from Laguerre *
C             polynomials.                                               *
C                                                                        *
C*************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=3)
      DIMENSION AINT(0:2*N-1),X(N),W(N),WK(N*N+4*N+1),IWK(N)
      DIMENSION WKD(N+1)
C
C     initialize; switch data if desired
C
      DATA AINT/1.D0,1.D0,2.D0,6.D0,24.D0,120.D0/
C
C     Output of test example
C
      WRITE (*,2000)
      WRITE (*,2010) (I,AINT(I),I=0,2*N-1)
C
C
      CALL ORTOGP (N,AINT,X,W,IERR,WK,WKD,IWK)
C
C
      IF (IERR.NE.0) THEN
C
C     Output of results
C
         WRITE (*,*)
     +    'ABANDONED COMPUTATIONS DUE TO ILL CONDITIONED MATRIX!'
      ELSE
         WRITE (*,2020)
         WRITE (*,2030) (I,X(I),W(I),I=1,N)
         WRITE (*,2040)
      ENDIF
      STOP
C
C
 2000 FORMAT (1X, 'C[ ', T74, ']*', /,
     +        1X, 'C[ TEST EXAMPLE:' ,T74, ']*', /,
     +        1X, 'C[ ', 13('='), T74, ']*', /,
     +        1X, 'C[', T74, ']*', /,
     +        1X, 'C[ GIVEN INTEGRAL VALUES:', 12X, 'I', 8X,
     +            'AINT(I)', T74, ']*', /,
     +        1X, 'C[ ', 30X, 25('-'), T74, ']*')
 2010 FORMAT (1X, 'C[ ', 30X, I5, '  I ', F15.5, T74, ']*')
 2020 FORMAT (1X, 'C[', T74, ']*', /,
     +        1X, 'C[ RESULTS:', T74, ']*', /,
     +        1X, 'C[ ', 4X, 'I', 8X, 'NODES X(I)', 8X,
     +            'WEIGHTS (I)', T74, ']*', /,
     +        1X, 'C[ ', 44('-'), T74, ']*')
 2030 FORMAT (1X, 'C[ ', I5,'  I ', F15.10, '  I ', F15.10, T74, ']*')
 2040 FORMAT (1X, 'C[', T74, ']*', /,
     +        1X, 'C[ ', 'NO ERROR', T74, ']*', /,
     +        1X, 'C[ ', T74, ']*')
      END
