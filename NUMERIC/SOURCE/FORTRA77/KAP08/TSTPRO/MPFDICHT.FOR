      PROGRAM TESFDI
C
C****************************************************************
C                                                               *
C  This program tests the subroutine  FDICHT  with functional   *
C  values of the sum of the complex basic functions             *
C                                                               *
C       2 * EXP(I*(-3)*OMEGA*X) + I * EXP(I*7*OMEGA*X)          *
C                                                               *
C  ( I : imaginary unit, OMEGA : 2*PI/P, P : period )           *
C  for different sets of equidistant nodes from the interval of *
C  periodicity  [0, P). It then compares the function values at *
C  the shifted nodes with the exact ones, which should differ   *
C  only the up to errors on the level of the machine constant.  * 
C                                                               *
C****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER M(5)
      DIMENSION FUNKRE(0:4095), FUNKIM(0:4095), P(5)
      DATA M / 16, 81, 263, 1101, 4096 /
      DATA P / 6.283185308D0, 7.0D0, 3.14159D0, 2.0D0, 10.0D0 /
      PI = 4.0D0*ATAN(1.0D0)
C      
      DO 50 I = 1, 5
C
C     Assign real and imaginary parts of the functional values
C
        FAKT = 2.0D0*PI/DBLE(M(I))
        DO 20 J = 0, M(I)-1
          XJ = DBLE(J)*FAKT
          FUNKRE(J) = 2.0D0*COS(3.0D0*XJ) - SIN(7.0D0*XJ)
          FUNKIM(J) = -2.0D0*SIN(3.0D0*XJ) + COS(7.0D0*XJ)
          FUNKIM(J)=0.0D0
  20    CONTINUE
        WEG = 0.5D0*P(I)/DBLE(M(I))
C
        CALL FDICHT (M(I), FUNKRE, FUNKIM, P(I), WEG)
C
C     Check the computed values
C
        FEHL = 0.0D0
        DO 30 J = 0, M(I)-1
          XJWEG = DBLE(J)*FAKT+2.0D0*PI/P(I)*WEG
          FRJ = 2.0D0*COS(3.0D0*XJWEG) - SIN(7.0D0*XJWEG)
          FIJ = -2.0D0*SIN(3.0D0*XJWEG) + COS(7.0D0*XJWEG)
          FIJ=0.0D0
          FEHLR = (FUNKRE(J) - FRJ)**2
          FEHLI = (FUNKIM(J) - FIJ)**2
          FEHL = MAX(FEHL,SQRT(FEHLR+FEHLI))
  30    CONTINUE
        WRITE(2,900) M(I), P(I), FEHL
  50  CONTINUE
      S T O P
 900  FORMAT(1X,'NUMBER = ',I4,'  PERIOD = ',F10.5,
     #       '  ERROR =  ',1PD12.5)
      END
