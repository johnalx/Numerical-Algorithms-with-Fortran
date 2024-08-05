      PROGRAM TEST
C                                       9.19.1992  (Dubois Guido)
C*****************************************************************
C                                                                *
C     Testprogram for the subroutine LAGUER.                     *
C     Compute the zeros of the polynomial                        *
C     PN(X) = A(N) * X**N +...+ A(1) * X + A(0)  with real roots *
C     using the method of Laguerre.                              *
C                                                                *
C     Test polynomial:                                           *
C     PN(X) = X**4 - 10 X**3 + 35 X**2 - 50 X + 24               *
C     Test results:                                              *
C                                                                *
C[    EXAMPLE OF SUBROUTINE LAGUER                              ]*
C[    ============================                              ]*
C[                                                              ]*
C[    COEFFICIENTS:                                             ]*
C[                                                              ]*
C[     I     A(I))                                              ]*
C[    --------------                                            ]*
C[     0     24.000                                             ]*
C[     1    -50.000                                             ]*
C[     2     35.000                                             ]*
C[     3    -10.000                                             ]*
C[     4      1.000                                             ]*
C[                                                              ]*
C[    ZEROS:                                                    ]*
C[                                                              ]*
C[     I     XI(I)    NITER                                     ]*
C[    ---------------------                                     ]*
C[     1      1.000     4                                       ]*
C[     2      2.000     4                                       ]*
C[     3      3.000     0                                       ]*
C[     4      4.000     0                                       ]*
C                                                                *
C*****************************************************************
C
      PARAMETER(N=4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(0:N),XI(1:N),WORK(0:N)
      INTEGER NITER(1:N)
C
C     Initialize; enter other polynomial if desired
C
      DATA A /24.,-50.,35.,-10.,1./
      DATA ABSERR /0./
      DATA RELERR /1.D-10/
      DATA IMAX /100/
C
C     Output of test coefficients
C
      WRITE(*,900)
      WRITE(*,910) (I,A(I),I=0,N)
      CALL LAGUER(A,N,ABSERR,RELERR,IMAX,XI,NITER,IANZ,WORK,IERR)
C
C     Output of test results
C
      IF(IERR .EQ. 1)THEN
         WRITE(*,920)
         WRITE(*,930) (I,XI(I),NITER(I),I=1,N)
      ELSE
         WRITE(*,940) IERR
      END IF
      STOP
  900 FORMAT(1X,'C[',4X,'EXAMPLE OF SUBROUTINE LAGUER',T66,']*',/,
     +       1X,'C[',4X,28('='),T66,']*',/,
     +       1X,'C[',T66,']*',/,
     +       1X,'C[',4X,'COEFFICIENTS:',T66,']*',/,
     +       1X,'C[',T66,']*',/,
     +       1X,'C[',4X,' I     A(I))',T66,']*',/,
     +       1X,'C[',4X,14('-'),T66,']*')
  910 FORMAT(1X,'C[',4X,I2,1X,F10.3,T66,']*')
  920 FORMAT(1X,'C[',T66,']*',/,
     +       1X,'C[',4X,'ZEROS:',T66,']*',/,
     +       1X,'C[',T66,']*',/,
     +       1X,'C[',4X,' I     XI(I)    NITER',T66,']*',/,
     +       1X,'C[',4X,21('-'),T66,']*')
  930 FORMAT(1X,'C[',4X,I2,1X,F10.3,3X,I3,T66,']*')
  940 FORMAT(1X,'C[',T66,']*',/,
     +       1X,'C[',4X,'ERROR WHEN EXECUTING SUBROUTINE LAGUER:',
     +       I2,T66,']*')
      END
