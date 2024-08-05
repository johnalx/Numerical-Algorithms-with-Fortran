      PROGRAM TEST
C                                       4.21.1992  (Dubois Guido)
C*****************************************************************
C                                                                *
C     Testprogram for the subroutine BAUPOL.                     *
C     We compute all real and complex roots of a polynomial      *
C     P(Z) = A(N) * Z**N +...+ A(1) * Z + A(0) with complex      *
C     coefficients.                                              *
C                                                                *
C     With the given polynomial                                  *
C     P(Z) = (Z-3*I-1) * (Z-2*I+1)**2 * (Z+4*I-2)                *
C     the test results are as follows:                           *
C                                                                *
C[    EXAMPLE OF SUBROUTINE BAUPOL                              ]*
C[    ============================                              ]*
C[                                                              ]*
C[    COEFFICIENTS:                                             ]*
C[                                                              ]*
C[     I   RE(A(I))   IM(A(I))                                  ]*
C[    -------------------------                                 ]*
C[     0  1.000E+00  0.000E+00                                  ]*
C[     1 -1.000E+00 -3.000E+00                                  ]*
C[     2  9.000E+00  1.200E+01                                  ]*
C[     3  4.900E+01 -4.300E+01                                  ]*
C[     4 -3.400E+01 -6.200E+01                                  ]*
C[                                                              ]*
C[    ZEROS:                                                    ]*
C[                                                              ]*
C[     I   RE(Z(I))   IM(Z(I))                                  ]*
C[    -------------------------                                 ]*
C[     1 -2.000E-01 -4.000E-01                                  ]*
C[     2  1.000E-01 -3.000E-01                                  ]*
C[     3  1.000E-01  2.000E-01                                  ]*
C[     4 -2.000E-01 -4.000E-01                                  ]*
C[                                                              ]*
C[       43 ITERATIONS                                          ]*
C                                                                *
C*****************************************************************
C
      PARAMETER(N=4)
      DOUBLE PRECISION KOEFRE(0:100),KOEFIM(0:100),
     +                 WURZRE(1:100),WURZIM(1:100)
C
C     Initialize; other test problems can be inserted here for 
C     differing N
C
      DATA KOEFRE/1.,-1., 9., 49.,-34.,96*0./
      DATA KOEFIM/0.,-3.,12.,-43.,-62.,96*0./
C
C     Output of Test coefficients
C
      WRITE(*,100)
      WRITE(*,110) (I,KOEFRE(I),KOEFIM(I),I=0,N)
      CALL BAUPOL(KOEFRE,KOEFIM,N,.TRUE.,WURZRE,WURZIM,ISHRIT)
C
C     Output of Test results
C
      WRITE(*,105)
      WRITE(*,110) (I,WURZRE(I),WURZIM(I),I=1,N)
      WRITE(*,120) ISHRIT
      STOP
  100 FORMAT(1X,'C[',4X,'EXAMPLE OF SUBROUTINE BAUPOL',T66,']*',/,
     +       1X,'C[',4X,28('='),T66,']*',/,
     +       1X,'C[',T66,']*',/,
     +       1X,'C[',4X,'COEFFICIENTS:',T66,']*',/,
     +       1X,'C[',T66,']*',/,
     +       1X,'C[',4X,' I   RE(A(I))   IM(A(I))',T66,']*',/,
     +       1X,'C[',4X,25('-'),T66,']*')
  110 FORMAT(1X,'C[',4X,I2,1X,1PE10.3,1X,E10.3,T66,']*')
  105 FORMAT(1X,'C[',T66,']*',/,
     +       1X,'C[',4X,'ZEROS:',T66,']*',/,
     +       1X,'C[',T66,']*',/,
     +       1X,'C[',4X,' I   RE(Z(I))   IM(Z(I))',T66,']*',/,
     +       1X,'C[',4X,25('-'),T66,']*')
  120 FORMAT(1X,'C[',T66,']*',/,
     +       1X,'C[',4X,I5,' ITERATIONS',T66,']*')
      END
